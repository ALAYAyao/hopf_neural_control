#include <px4_platform_common/module.h>
#include <px4_platform_common/module_params.h>
#include <px4_platform_common/log.h>
#include <px4_platform_common/tasks.h>

#include <lib/perf/perf_counter.h>

#include <drivers/drv_hrt.h>

#include <uORB/Publication.hpp>
#include <uORB/topics/auv_hopf_status.h>
#include <uORB/topics/debug_key_value.h>
#include <string.h>

#include <math.h>

#include "auv_common.hpp"
#include "auv_model.hpp"
#include "hopf_observer.hpp"
#include "wave_disturbance.hpp"
#include "mpc_optim.hpp"
#include "actuator_saturation.hpp"
#include "telem_link.hpp"

/**
 * @file auv_hopf_main.cpp
 *
 * @brief AUV Hopf + Tube-MPC 主模块
 *
 * 当前版本定位：
 * 1. 先把学术代码主链在 PX4 模块内部完整跑通
 * 2. 使用内部状态 x / z 做闭环仿真
 * 3. 明确拆出第一层 central-path OCP 与第二层扰动增强预测层
 *
 * 主要数据流（当前两层过渡版）：
 *   mission_ref_hor = build_mission_reference_horizon(t)
 *   d_true          = wave_disturbance.step(x, t)
 *   d_hat           = hopf_observer.update(x, tau_prev, dt)
 *
 *   v_central       = central_mpc.solve(z, mission_ref_hor, d_zero)
 *   z_hor, v_hor    = central_mpc.export_last_solution_horizon(...)
 *
 *   tau_anc         = ancillary_mpc.solve(x, z_hor, v_hor, d_hat_hor)
 *   tau_applied     = actuator_saturation.apply(tau_anc)
 *
 *   x_next          = model.step(x, tau_applied, d_true, dt)
 *   z_next          = model.step(z, v_central, 0, dt)
 *
 * 注意：
 * - 手工圆轨迹不再直接扮演第一层 nominal trajectory；
 * - 它现在只是 mission reference，真正的第一层由 central_mpc 在线生成；
 * - 第二层现在优先使用第一层最近一次求解导出的整条 z* / v* 时域；
 * - 若导出失败，才退回到近似前向滚动兜底；
 * - 执行层现在显式使用 actuator_saturation，thruster_allocation 不再承担控制器语义。
 * - 参数初始化改为通过 profile bundle 一次性装载，正式区分 paper_reproduction 与 engineering_default。
 */

class AuvHopf :
    public ModuleBase,
    public ModuleParams
{
public:
    static Descriptor desc;

    AuvHopf();
    ~AuvHopf() override;

    static int trampoline(int argc, char *argv[]);
    static int task_spawn(int argc, char *argv[]);
    static AuvHopf *instantiate(int argc, char *argv[]);
    static int custom_command(int argc, char *argv[]);
    static int print_usage(const char *reason = nullptr);

    bool init();
    int print_status() override;
    void run() override;

private:
    // --------------------------------------------------
    // 核心算法对象
    // --------------------------------------------------
    AuvModel _model;
    HopfObserver _hopf;
    WaveDisturbance _wave;
    MpcOptim _central_mpc;
    MpcOptim _anc_mpc;
    ActuatorSaturation _actuator_sat;
    TelemLink _telem{"/dev/ttyS1"};

    //----------------------------------------------
    // 性能计数器
    //----------------------------------------------
    perf_counter_t _perf_central_solve{perf_alloc(PC_ELAPSED, "auv_hopf_central_mpc")};
    perf_counter_t _perf_anc_solve{perf_alloc(PC_ELAPSED, "auv_hopf_anc_mpc")};
    perf_counter_t _perf_total_run{perf_alloc(PC_ELAPSED, "auv_hopf_total")};

    // --------------------------------------------------
    // 状态变量
    // x : 真实状态 actual state
    // z : central-path 当前名义状态
    // --------------------------------------------------
    auv::StateVec _x{};
    auv::StateVec _z{};

    // mission-level 当前参考点（上层任务参考，不直接等于第一层）
    auv::StateVec _xref{};

    // 当前真实扰动 / 学习估计 / 控制量
    auv::InputVec _d_true{};
    auv::InputVec _d_hat{};
    auv::InputVec _v_central{};
    auv::InputVec _tau_anc{};
    auv::InputVec _tau_applied{};
    auv::InputVec _tau_prev{};
    auv::InputVec _zero_input{};

    // 调度时间
    hrt_abstime _last_run{0};
    float _sim_time{0.0f};

    // 时域缓存
    float _mission_hor[auv::NX * auv::NHOR_MAX]{};
    float _zref_hor[auv::NX * auv::NHOR_MAX]{};
    float _vref_hor[auv::NU * auv::NHOR_MAX]{};
    float _dhat_hor[auv::NP * auv::NHOR_MAX]{};
    float _dzero_hor[auv::NP * auv::NHOR_MAX]{};
    int _xref_nhor{0};

    // 分层约束缓存：
    // nominal central layer / ancillary layer / physical layer
    auv::MpcStateBoundsParams _central_state_bounds{};
    auv::MpcStateBoundsParams _anc_state_bounds{};

    // 参数口径与任务配置
    static constexpr auv::ParameterProfile kParameterProfile = auv::ParameterProfile::EngineeringDefault;
    const char *_parameter_profile_name{auv::parameter_profile_name(kParameterProfile)};
    auv::MissionReferenceParams _mission_ref_params{};
    auv::ExperimentConfig _experiment_cfg{};
    auv::PhysicalActuatorLimits _physical_limits{};

    // 理论-实现对应的在线指标
    float _dhat_prev_hor[auv::NP * auv::NHOR_MAX]{};
    float _disturbance_residual_norm{0.0f};
    float _preview_update_norm{0.0f};
    float _tube_tracking_error_norm{0.0f};
    float _tube_value_surrogate{0.0f};
    float _tube_value_delta{0.0f};
    float _central_state_violation_norm{0.0f};
    float _anc_state_violation_norm{0.0f};
    float _central_input_violation_norm{0.0f};
    float _anc_input_violation_norm{0.0f};
    float _physical_clip_norm{0.0f};
    bool _central_export_ok{false};
    bool _central_bridge_ok{true};
    bool _anc_bridge_ok{true};

    // 调试输出
    uORB::Publication<auv_hopf_status_s> _status_pub{ORB_ID(auv_hopf_status)};
    uORB::Publication<debug_key_value_s> _debug_key_value_pub{ORB_ID(debug_key_value)};

    // 内部帮助函数
    void reset_states();
    void build_mission_reference(float t);
    void build_mission_reference_horizon();
    void build_disturbance_preview_horizon();
    void build_zero_disturbance_horizon();
    void build_central_rollout_horizon();
    void clamp_state_horizon_to_bounds(float *x_hor,
                                       int nhor,
                                       const auv::MpcStateBoundsParams &bounds);
    void publish_debug(float dt);
    void publish_scalar_debug(const char *key, float value);

    static float sanitize_dt(float dt);
    static float clampf(float v, float lo, float hi);
    static float wrap_pi(float a);
    static float input_norm(const auv::InputVec &v);
    static float tube_error_norm(const auv::StateVec &x, const auv::StateVec &z);
    static float state_bound_violation_norm(const auv::StateVec &x,
                                            const auv::MpcStateBoundsParams &bounds);
    static float input_bound_violation_norm(const auv::InputVec &u,
                                            const auv::MpcBoundsParams &bounds);
    static float quadratic_tracking_value(const auv::StateVec &x,
                                          const auv::StateVec &z,
                                          const auv::MpcWeightParams &weights);
};

AuvHopf::AuvHopf() :
    ModuleParams(nullptr)
{
    reset_states();
}

AuvHopf::~AuvHopf()
{
    perf_free(_perf_central_solve);
    perf_free(_perf_anc_solve);
    perf_free(_perf_total_run);
}

float AuvHopf::clampf(float v, float lo, float hi)
{
    if (v < lo) {
        return lo;
    }

    if (v > hi) {
        return hi;
    }

    return v;
}

float AuvHopf::wrap_pi(float a)
{
    while (a > M_PI_F) {
        a -= 2.0f * M_PI_F;
    }

    while (a < -M_PI_F) {
        a += 2.0f * M_PI_F;
    }

    return a;
}

float AuvHopf::input_norm(const auv::InputVec &v)
{
    return sqrtf(v(0) * v(0) + v(1) * v(1) + v(2) * v(2));
}

float AuvHopf::tube_error_norm(const auv::StateVec &x, const auv::StateVec &z)
{
    const float e0 = x(0) - z(0);
    const float e1 = x(1) - z(1);
    const float e2 = wrap_pi(x(2) - z(2));
    const float e3 = x(3) - z(3);
    const float e4 = x(4) - z(4);
    const float e5 = x(5) - z(5);

    return sqrtf(e0 * e0 + e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4 + e5 * e5);
}

float AuvHopf::state_bound_violation_norm(const auv::StateVec &x,
                                          const auv::MpcStateBoundsParams &bounds)
{
    if (!bounds.enable_soft_bounds) {
        return 0.0f;
    }

    float accum = 0.0f;

    for (int i = 0; i < auv::NX; ++i) {
        float viol = 0.0f;

        if (x(i) < bounds.xmin[i]) {
            viol = bounds.xmin[i] - x(i);
        } else if (x(i) > bounds.xmax[i]) {
            viol = x(i) - bounds.xmax[i];
        }

        accum += viol * viol;
    }

    return sqrtf(accum);
}

float AuvHopf::input_bound_violation_norm(const auv::InputVec &u,
                                          const auv::MpcBoundsParams &bounds)
{
    float accum = 0.0f;

    for (int i = 0; i < auv::NU; ++i) {
        float viol = 0.0f;

        if (u(i) < bounds.umin[i]) {
            viol = bounds.umin[i] - u(i);
        } else if (u(i) > bounds.umax[i]) {
            viol = u(i) - bounds.umax[i];
        }

        accum += viol * viol;
    }

    return sqrtf(accum);
}

float AuvHopf::quadratic_tracking_value(const auv::StateVec &x,
                                        const auv::StateVec &z,
                                        const auv::MpcWeightParams &weights)
{
    float value = 0.0f;

    for (int i = 0; i < auv::NX; ++i) {
        float e = x(i) - z(i);
        if (i == 2) {
            e = wrap_pi(e);
        }

        value += 0.5f * weights.Q[i] * e * e;
    }

    return value;
}

float AuvHopf::sanitize_dt(float dt)
{
    if (!PX4_ISFINITE(dt) || dt <= 0.0f) {
        return 0.02f;
    }

    if (dt > 0.2f) {
        return 0.2f;
    }

    return dt;
}

void AuvHopf::reset_states()
{
    _x.zero();
    _z.zero();
    _xref.zero();

    _d_true.zero();
    _d_hat.zero();
    _v_central.zero();
    _tau_anc.zero();
    _tau_applied.zero();
    _tau_prev.zero();
    _zero_input.zero();

    _sim_time = 0.0f;
    _last_run = 0;
    _xref_nhor = 0;
    _disturbance_residual_norm = 0.0f;
    _preview_update_norm = 0.0f;
    _tube_tracking_error_norm = 0.0f;
    _tube_value_surrogate = 0.0f;
    _tube_value_delta = 0.0f;
    _central_state_violation_norm = 0.0f;
    _anc_state_violation_norm = 0.0f;
    _central_input_violation_norm = 0.0f;
    _anc_input_violation_norm = 0.0f;
    _physical_clip_norm = 0.0f;
    _central_export_ok = false;
    _central_bridge_ok = true;
    _anc_bridge_ok = true;

    _central_state_bounds = auv::default_central_state_bounds_params();
    _anc_state_bounds = auv::default_ancillary_state_bounds_params();

    for (int i = 0; i < auv::NX * auv::NHOR_MAX; ++i) {
        _mission_hor[i] = 0.0f;
        _zref_hor[i] = 0.0f;
    }

    for (int i = 0; i < auv::NU * auv::NHOR_MAX; ++i) {
        _vref_hor[i] = 0.0f;
    }

    for (int i = 0; i < auv::NP * auv::NHOR_MAX; ++i) {
        _dhat_hor[i] = 0.0f;
        _dhat_prev_hor[i] = 0.0f;
        _dzero_hor[i] = 0.0f;
    }

    // 初始真实状态
    _x(0) = 0.0f;
    _x(1) = 0.0f;
    _x(2) = 0.0f;
    _x(3) = 0.0f;
    _x(4) = 0.0f;
    _x(5) = 0.0f;

    // 初始名义状态
    _z = _x;

    _wave.reset();
    _hopf.reset();
}

/**
 * @brief 构造当前 mission reference 点
 *
 * 这里仍以平面圆轨迹为例，但它现在只承担“任务参考”角色：
 *   Xref = R cos(wt)
 *   Yref = R sin(wt)
 *   psiref = atan2(Ydot, Xdot)
 */
void AuvHopf::build_mission_reference(float t)
{
    const float R = _mission_ref_params.circle_radius;
    const float w = _mission_ref_params.circle_omega;

    const float xr = R * cosf(w * t);
    const float yr = R * sinf(w * t);

    const float xdot_r = -R * w * sinf(w * t);
    const float ydot_r =  R * w * cosf(w * t);

    const float psi_r = atan2f(ydot_r, xdot_r);
    const float u_r = R * w;
    const float v_r = 0.0f;
    const float r_r = w;

    _xref.zero();
    _xref(0) = xr;
    _xref(1) = yr;
    _xref(2) = psi_r;
    _xref(3) = u_r;
    _xref(4) = v_r;
    _xref(5) = r_r;
}

/**
 * @brief 构造 mission reference horizon
 *
 * 当前手工圆轨迹只作为第一层 nominal OCP 的外部目标，
 * 不再直接被解释成论文里的 central path。
 */
void AuvHopf::build_mission_reference_horizon()
{
    _xref_nhor = (_central_mpc.horizon_length() > 0) ? _central_mpc.horizon_length() : 20;
    const float dt = (_central_mpc.sample_time() > 0.0f) ? _central_mpc.sample_time() : 0.02f;

    for (int k = 0; k < auv::NHOR_MAX; ++k) {
        for (int i = 0; i < auv::NX; ++i) {
            _mission_hor[k * auv::NX + i] = 0.0f;
        }
    }

    for (int k = 0; k < _xref_nhor && k < auv::NHOR_MAX; ++k) {
        const float tk = _sim_time + static_cast<float>(k) * dt;

        const float R = _mission_ref_params.circle_radius;
        const float w = _mission_ref_params.circle_omega;

        const float xr = R * cosf(w * tk);
        const float yr = R * sinf(w * tk);

        const float xdot_r = -R * w * sinf(w * tk);
        const float ydot_r =  R * w * cosf(w * tk);

        const float psi_r = atan2f(ydot_r, xdot_r);
        const float u_r = R * w;
        const float v_r = 0.0f;
        const float r_r = w;

        _mission_hor[k * auv::NX + 0] = xr;
        _mission_hor[k * auv::NX + 1] = yr;
        _mission_hor[k * auv::NX + 2] = psi_r;
        _mission_hor[k * auv::NX + 3] = u_r;
        _mission_hor[k * auv::NX + 4] = v_r;
        _mission_hor[k * auv::NX + 5] = r_r;
    }

    clamp_state_horizon_to_bounds(_mission_hor, _xref_nhor, _central_state_bounds);
}

void AuvHopf::build_disturbance_preview_horizon()
{
    const int nh = (_xref_nhor > 1) ? _xref_nhor : 2;

    for (int k = 0; k < auv::NHOR_MAX; ++k) {
        for (int i = 0; i < auv::NP; ++i) {
            _dhat_hor[k * auv::NP + i] = 0.0f;
        }
    }

    // 当前仍采用 zero-order hold preview：
    // 先保证第二层预测模型真的能“看到” Hopf 估计。
    for (int k = 0; k < nh && k < auv::NHOR_MAX; ++k) {
        _dhat_hor[k * auv::NP + 0] = _d_hat(0);
        _dhat_hor[k * auv::NP + 1] = _d_hat(1);
        _dhat_hor[k * auv::NP + 2] = _d_hat(2);
    }

    for (int k = nh; k < auv::NHOR_MAX; ++k) {
        _dhat_hor[k * auv::NP + 0] = _dhat_hor[(nh - 1) * auv::NP + 0];
        _dhat_hor[k * auv::NP + 1] = _dhat_hor[(nh - 1) * auv::NP + 1];
        _dhat_hor[k * auv::NP + 2] = _dhat_hor[(nh - 1) * auv::NP + 2];
    }

    _preview_update_norm = 0.0f;
    for (int k = 0; k < nh && k < auv::NHOR_MAX; ++k) {
        const float du = _dhat_hor[k * auv::NP + 0] - _dhat_prev_hor[k * auv::NP + 0];
        const float dv = _dhat_hor[k * auv::NP + 1] - _dhat_prev_hor[k * auv::NP + 1];
        const float dr = _dhat_hor[k * auv::NP + 2] - _dhat_prev_hor[k * auv::NP + 2];
        const float delta = sqrtf(du * du + dv * dv + dr * dr);
        if (delta > _preview_update_norm) {
            _preview_update_norm = delta;
        }
    }

    for (int i = 0; i < auv::NP * auv::NHOR_MAX; ++i) {
        _dhat_prev_hor[i] = _dhat_hor[i];
    }
}

void AuvHopf::build_zero_disturbance_horizon()
{
    for (int k = 0; k < auv::NHOR_MAX; ++k) {
        for (int i = 0; i < auv::NP; ++i) {
            _dzero_hor[k * auv::NP + i] = 0.0f;
        }
    }
}

void AuvHopf::clamp_state_horizon_to_bounds(float *x_hor,
                                            int nhor,
                                            const auv::MpcStateBoundsParams &bounds)
{
    if (x_hor == nullptr || !bounds.enable_soft_bounds) {
        return;
    }

    const int nh = (nhor < 2) ? 2 : ((nhor > auv::NHOR_MAX) ? auv::NHOR_MAX : nhor);

    for (int k = 0; k < nh; ++k) {
        for (int i = 0; i < auv::NX; ++i) {
            if (bounds.penalty[i] <= 0.0f) {
                continue;
            }

            x_hor[k * auv::NX + i] = clampf(x_hor[k * auv::NX + i],
                                            bounds.xmin[i],
                                            bounds.xmax[i]);
        }
    }
}

/**
 * @brief 获取第一层 recent solve 对应的 z* / v* horizon
 *
 * 新版优先直接读取 central_mpc 最近一次求解得到的整条最优控制时域，
 * 再由求解器侧缓存出的状态/控制时域作为第二层参考。
 *
 * 只有在导出失败时，才退回到“当前 z + 当前 v_central 的前向滚动”近似。
 */
void AuvHopf::build_central_rollout_horizon()
{
    const int nh = (_xref_nhor > 1) ? _xref_nhor : 2;

    for (int k = 0; k < auv::NHOR_MAX; ++k) {
        for (int i = 0; i < auv::NX; ++i) {
            _zref_hor[k * auv::NX + i] = 0.0f;
        }

        for (int i = 0; i < auv::NU; ++i) {
            _vref_hor[k * auv::NU + i] = 0.0f;
        }
    }

    _central_export_ok = _central_mpc.export_last_solution_horizon(_zref_hor, _vref_hor, nh);
    if (_central_export_ok) {
        return;
    }

    const float dt = (_central_mpc.sample_time() > 0.0f) ? _central_mpc.sample_time() : 0.02f;

    auv::StateVec z_roll = _z;
    auv::InputVec d_zero{};
    d_zero.zero();

    for (int k = 0; k < nh && k < auv::NHOR_MAX; ++k) {
        for (int i = 0; i < auv::NX; ++i) {
            _zref_hor[k * auv::NX + i] = z_roll(i);
        }

        for (int i = 0; i < auv::NU; ++i) {
            _vref_hor[k * auv::NU + i] = _v_central(i);
        }

        z_roll = _model.step(z_roll, _v_central, d_zero, dt);
    }

    for (int k = nh; k < auv::NHOR_MAX; ++k) {
        for (int i = 0; i < auv::NX; ++i) {
            _zref_hor[k * auv::NX + i] = _zref_hor[(nh - 1) * auv::NX + i];
        }

        for (int i = 0; i < auv::NU; ++i) {
            _vref_hor[k * auv::NU + i] = _vref_hor[(nh - 1) * auv::NU + i];
        }
    }

    clamp_state_horizon_to_bounds(_zref_hor, nh, _central_state_bounds);
}

bool AuvHopf::init()
{
    const auto bundle = auv::make_runtime_profile_bundle(kParameterProfile);

    _parameter_profile_name = bundle.profile_name;
    _mission_ref_params = bundle.mission_ref;
    _experiment_cfg = bundle.experiment;
    _physical_limits = bundle.physical_limits;

    _model.set_plant_params(bundle.plant);

    _hopf.set_plant_params(bundle.plant);
    _hopf.set_hopf_config(bundle.hopf);

    _wave.set_wave_env(bundle.wave_env);
    _wave.set_model_params(bundle.wave_model);

    _actuator_sat.set_limits(bundle.physical_limits);

    _central_state_bounds = bundle.central_state_bounds;
    _anc_state_bounds = bundle.ancillary_state_bounds;

    _central_mpc.set_plant_params(bundle.plant);
    _central_mpc.set_weight_params(bundle.central_weights);
    _central_mpc.set_bounds_params(bundle.central_bounds);
    _central_mpc.set_state_bounds_params(bundle.central_state_bounds);
    _central_mpc.set_mpc_config(bundle.central_config);
    _central_mpc.initialize(bundle.central_config.dt, bundle.central_config.Thor);

    _anc_mpc.set_plant_params(bundle.plant);
    _anc_mpc.set_weight_params(bundle.ancillary_weights);
    _anc_mpc.set_bounds_params(bundle.ancillary_bounds);
    _anc_mpc.set_state_bounds_params(bundle.ancillary_state_bounds);
    _anc_mpc.set_mpc_config(bundle.ancillary_config);
    _anc_mpc.initialize(bundle.ancillary_config.dt, bundle.ancillary_config.Thor);

    reset_states();

    _telem.init();

    PX4_INFO("auv_hopf started, profile=%s", _parameter_profile_name);
    return true;
}

void AuvHopf::run()
{
    if (!init()) {
        PX4_ERR("init failed");
        return;
    }

    float dt = (_central_mpc.sample_time() > 0.0f) ? _central_mpc.sample_time() : 0.02f;
    dt = sanitize_dt(dt);
    const hrt_abstime dt_us = static_cast<hrt_abstime>(dt * 1e6f);

    hrt_abstime next_run_time = hrt_absolute_time();

    while (!should_exit()) {

        next_run_time += dt_us;

        perf_begin(_perf_total_run);
        _sim_time += dt;

        // --------------------------------------------------
        // 1) 构造 mission reference 与外扰 / 学习器输出
        build_mission_reference(_sim_time);
        build_mission_reference_horizon();

        _d_true = _wave.step(_x, _sim_time);
        _d_hat  = _hopf.update(_x, _tau_prev, dt);
        _disturbance_residual_norm = input_norm(_d_true - _d_hat);

        build_disturbance_preview_horizon();
        build_zero_disturbance_horizon();

        // --------------------------------------------------
        // 2) 第一层：central-path OCP
        //
        // 目标：
        //   z, v 由 nominal OCP 追 mission reference 生成；
        //   这里显式不注入 Hopf preview，保持第一层干净。
        _central_mpc.set_reference_horizon_control_preview(
            _z,
            _mission_hor,
            nullptr,
            _dzero_hor,
            _xref_nhor);

        perf_begin(_perf_central_solve);
        _v_central = _central_mpc.solve();
        perf_end(_perf_central_solve);
        _central_bridge_ok = _central_mpc.debug().bridge_owner_match;
        _central_input_violation_norm = input_bound_violation_norm(_v_central,
                                                                   _central_mpc.bounds_params());

        // --------------------------------------------------
        // 3) 用 central 输出构造 z* / v* horizon，供第二层使用
        build_central_rollout_horizon();

        // --------------------------------------------------
        // 4) 第二层：Hopf-augmented ancillary MPC
        //
        // 当前第二层使用：
        //   xref = z* horizon
        //   uref = v* horizon
        //   p    = d_hat preview
        _anc_mpc.set_reference_horizon_control_preview(
            _x,
            _zref_hor,
            _vref_hor,
            _dhat_hor,
            _xref_nhor);

        perf_begin(_perf_anc_solve);
        _tau_anc = _anc_mpc.solve();
        perf_end(_perf_anc_solve);
        _anc_bridge_ok = _anc_mpc.debug().bridge_owner_match;
        _anc_input_violation_norm = input_bound_violation_norm(_tau_anc,
                                                               _anc_mpc.bounds_params());

        // --------------------------------------------------
        // 5) 执行层（当前过渡版）
        //
        // 先借用 thruster_allocation 只做物理限幅，不再让它承担
        // “第二层 controller”职责。
        auv::InputVec d_zero{};
        d_zero.zero();
        _tau_applied = _actuator_sat.apply(_tau_anc);
        _physical_clip_norm = input_norm(_actuator_sat.debug().tau_clip);

        // --------------------------------------------------
        // 6) 推进真实状态与 central 状态
        _x = _model.step(_x, _tau_applied, _d_true, dt);
        _z = _model.step(_z, _v_central, d_zero, dt);

        // --------------------------------------------------
        // 7) 记录上一拍实际执行控制
        _tau_prev = _tau_applied;
        _tube_tracking_error_norm = tube_error_norm(_x, _z);
        const float tube_value_prev = _tube_value_surrogate;
        _tube_value_surrogate = quadratic_tracking_value(_x, _z, _anc_mpc.weight_params());
        _tube_value_delta = _tube_value_surrogate - tube_value_prev;
        _central_state_violation_norm = state_bound_violation_norm(_z, _central_state_bounds);
        _anc_state_violation_norm = state_bound_violation_norm(_x, _anc_state_bounds);

        // --------------------------------------------------
        // 8) 发布调试信息
        publish_debug(dt);

        perf_end(_perf_total_run);

        hrt_abstime now = hrt_absolute_time();

        if (_sim_time >= _experiment_cfg.sim_time_limit) {
            PX4_INFO("时间到！已运行 %.1f 秒，算法自动安全退出。", (double)_experiment_cfg.sim_time_limit);
            request_stop();
        }

        if (now < next_run_time) {
            px4_usleep(next_run_time - now);
        } else {
            PX4_WARN("Deadline Missed! dt overrun.");
            next_run_time = now;
        }
    }
}

int AuvHopf::print_status()
{
    PX4_INFO("running");
    PX4_INFO("profile = %s", _parameter_profile_name);
    PX4_INFO("t = %.3f", (double)_sim_time);

    PX4_INFO("x = [%.3f %.3f %.3f %.3f %.3f %.3f]",
             (double)_x(0), (double)_x(1), (double)_x(2),
             (double)_x(3), (double)_x(4), (double)_x(5));

    PX4_INFO("z = [%.3f %.3f %.3f %.3f %.3f %.3f]",
             (double)_z(0), (double)_z(1), (double)_z(2),
             (double)_z(3), (double)_z(4), (double)_z(5));

    PX4_INFO("mission_ref = [%.3f %.3f %.3f %.3f %.3f %.3f]",
             (double)_xref(0), (double)_xref(1), (double)_xref(2),
             (double)_xref(3), (double)_xref(4), (double)_xref(5));

    PX4_INFO("x-z = [%.3f %.3f %.3f %.3f %.3f %.3f]",
             (double)(_x(0) - _z(0)),
             (double)(_x(1) - _z(1)),
             (double)(_x(2) - _z(2)),
             (double)(_x(3) - _z(3)),
             (double)(_x(4) - _z(4)),
             (double)(_x(5) - _z(5)));

    PX4_INFO("d_true = [%.3f %.3f %.3f]",
             (double)_d_true(0), (double)_d_true(1), (double)_d_true(2));

    PX4_INFO("d_hat = [%.3f %.3f %.3f]",
             (double)_d_hat(0), (double)_d_hat(1), (double)_d_hat(2));

    PX4_INFO("v_central = [%.3f %.3f %.3f]",
             (double)_v_central(0), (double)_v_central(1), (double)_v_central(2));

    PX4_INFO("tau_anc = [%.3f %.3f %.3f]",
             (double)_tau_anc(0), (double)_tau_anc(1), (double)_tau_anc(2));

    PX4_INFO("tau_applied = [%.3f %.3f %.3f]",
             (double)_tau_applied(0), (double)_tau_applied(1), (double)_tau_applied(2));

    PX4_INFO("metrics = residual %.3f, preview_delta %.3f, tube_error %.3f, tube_value %.3f, tube_dV %.3f",
             (double)_disturbance_residual_norm,
             (double)_preview_update_norm,
             (double)_tube_tracking_error_norm,
             (double)_tube_value_surrogate,
             (double)_tube_value_delta);

    PX4_INFO("violations = z_state %.3f, x_state %.3f, v_central %.3f, tau_anc %.3f, clip %.3f",
             (double)_central_state_violation_norm,
             (double)_anc_state_violation_norm,
             (double)_central_input_violation_norm,
             (double)_anc_input_violation_norm,
             (double)_physical_clip_norm);

    PX4_INFO("solver_guard = central_bridge %d, ancillary_bridge %d, central_export %d, central_stub %d, ancillary_stub %d",
             _central_bridge_ok ? 1 : 0,
             _anc_bridge_ok ? 1 : 0,
             _central_export_ok ? 1 : 0,
             _central_mpc.debug().rodas_stub_touch_count,
             _anc_mpc.debug().rodas_stub_touch_count);

    PX4_INFO("xref_nhor = %d", _xref_nhor);
    PX4_INFO("preview0 = [%.3f %.3f %.3f]",
             (double)_dhat_hor[0],
             (double)_dhat_hor[1],
             (double)_dhat_hor[2]);

    PX4_INFO("zref0 = [%.3f %.3f %.3f %.3f %.3f %.3f]",
             (double)_zref_hor[0], (double)_zref_hor[1], (double)_zref_hor[2],
             (double)_zref_hor[3], (double)_zref_hor[4], (double)_zref_hor[5]);

    return 0;
}

AuvHopf *AuvHopf::instantiate(int argc, char *argv[])
{
    return new AuvHopf();
}

// ==========================================================
// 线程入口跳板：系统在这里创建实例，并卡入死循环
// ==========================================================
int AuvHopf::trampoline(int argc, char *argv[])
{
    AuvHopf *instance = new AuvHopf();
    if (!instance) {
        PX4_ERR("alloc failed");
        return -1;
    }

    desc.object.store(instance);
    instance->run();
    delete instance;
    desc.object.store(nullptr);
    desc.task_id = -1;
    return 0;
}

int AuvHopf::task_spawn(int argc, char *argv[])
{
    int task_id = px4_task_spawn_cmd("auv_hopf",
                                     SCHED_DEFAULT,
                                     SCHED_PRIORITY_ATTITUDE_CONTROL,
                                     4096,
                                     (px4_main_t)&AuvHopf::trampoline,
                                     (char *const *)argv);

    if (task_id < 0) {
        return -errno;
    }

    desc.task_id = task_id;
    return 0;
}

int AuvHopf::custom_command(int argc, char *argv[])
{
    return print_usage("unknown command");
}

int AuvHopf::print_usage(const char *reason)
{
    if (reason) {
        PX4_WARN("%s", reason);
    }

    PRINT_MODULE_DESCRIPTION(
        R"DESCR_STR(
### Description
AUV Hopf + Tube-MPC composite control module.

Current version runs a two-layer transitional architecture internally:
wave disturbance -> Hopf learner -> central-path MPC -> Hopf-augmented ancillary MPC -> physical saturation -> plant update

Key property of this version:
- mission reference is no longer treated as the central path itself
- central state z is generated by a dedicated nominal OCP
- ancillary layer explicitly tracks z* / v* and uses disturbance preview
- actuator_saturation is now the dedicated physical layer object
- runtime parameters are loaded through a named profile bundle: engineering_default or paper_reproduction
- debug stream now exposes residual / preview-update / tube-value / bound-violation metrics for formal validation
)DESCR_STR");

    PRINT_MODULE_USAGE_NAME("auv_hopf", "template");
    PRINT_MODULE_USAGE_COMMAND("start");
    PRINT_MODULE_USAGE_COMMAND("status");
    PRINT_MODULE_USAGE_COMMAND("stop");
    return 0;
}

void AuvHopf::publish_debug(float dt)
{
    auv_hopf_status_s msg{};
    msg.timestamp = hrt_absolute_time();

    msg.sim_time = _sim_time;
    msg.dt = dt;

    for (int i = 0; i < auv::NX; ++i) {
        msg.x[i]         = _x(i);
        msg.z[i]         = _z(i);
        msg.xref[i]      = _xref(i);
        msg.error_x_z[i] = _x(i) - _z(i);
    }

    for (int i = 0; i < auv::NU; ++i) {
        // 这里先把 tau_nom 字段复用为第一层 central input，
        // 便于后续区分第一层 / 第二层控制输出。
        msg.d_true[i]      = _d_true(i);
        msg.d_hat[i]       = _d_hat(i);
        msg.tau_nom[i]     = _v_central(i);
        msg.tau_applied[i] = _tau_applied(i);
    }

    const auto &hopf_dbg = _hopf.debug();
    msg.hopf_surge[0] = hopf_dbg.surge.d_raw;
    msg.hopf_surge[1] = hopf_dbg.surge.d_pred_before;
    msg.hopf_surge[2] = hopf_dbg.surge.f_error;

    msg.hopf_sway[0] = hopf_dbg.sway.d_raw;
    msg.hopf_sway[1] = hopf_dbg.sway.d_pred_before;
    msg.hopf_sway[2] = hopf_dbg.sway.f_error;

    msg.hopf_yaw[0] = hopf_dbg.yaw.d_raw;
    msg.hopf_yaw[1] = hopf_dbg.yaw.d_pred_before;
    msg.hopf_yaw[2] = hopf_dbg.yaw.f_error;

    const auto &sat_dbg = _actuator_sat.debug();
    for (int i = 0; i < 3; ++i) {
        msg.tau_fb[i]        = 0.0f;
        msg.tau_ancillary[i] = _tau_anc(i);
        msg.tau_raw[i]       = sat_dbg.tau_cmd(i);
    }

    msg.xref_nhor = static_cast<uint8_t>(_xref_nhor);
    _status_pub.publish(msg);

    publish_scalar_debug("resid", _disturbance_residual_norm);
    publish_scalar_debug("prvupd", _preview_update_norm);
    publish_scalar_debug("tubee", _tube_tracking_error_norm);
    publish_scalar_debug("tubev", _tube_value_surrogate);
    publish_scalar_debug("tubedv", _tube_value_delta);
    publish_scalar_debug("zviol", _central_state_violation_norm);
    publish_scalar_debug("xviol", _anc_state_violation_norm);
    publish_scalar_debug("cvu", _central_input_violation_norm);
    publish_scalar_debug("avu", _anc_input_violation_norm);
    publish_scalar_debug("clipn", _physical_clip_norm);
    publish_scalar_debug("cexok", _central_export_ok ? 1.0f : 0.0f);
    publish_scalar_debug("cbrdg", _central_bridge_ok ? 1.0f : 0.0f);
    publish_scalar_debug("abrdg", _anc_bridge_ok ? 1.0f : 0.0f);
    publish_scalar_debug("clipu", sat_dbg.tau_clip(0));
    publish_scalar_debug("clipv", sat_dbg.tau_clip(1));
    publish_scalar_debug("clipr", sat_dbg.tau_clip(2));
}

void AuvHopf::publish_scalar_debug(const char *key, float value)
{
    debug_key_value_s msg{};
    msg.timestamp = hrt_absolute_time();
    msg.value = value;
    strncpy(msg.key, key, sizeof(msg.key) - 1);
    msg.key[sizeof(msg.key) - 1] = '\0';
    _debug_key_value_pub.publish(msg);
}

ModuleBase::Descriptor AuvHopf::desc{
    AuvHopf::task_spawn, AuvHopf::custom_command, AuvHopf::print_usage
};

extern "C" __EXPORT int auv_hopf_main(int argc, char *argv[])
{
    return ModuleBase::main(AuvHopf::desc, argc, argv);
}
