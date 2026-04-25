#include <px4_platform_common/module.h>
#include <px4_platform_common/module_params.h>
#include <px4_platform_common/log.h>
#include <px4_platform_common/tasks.h>

#include <lib/perf/perf_counter.h>

#include <drivers/drv_hrt.h>

#include <uORB/Publication.hpp>
#include <uORB/topics/auv_hopf_status.h>
#include <string.h>

#include <math.h>
#include <errno.h>

#include "auv_common.hpp"
#include "auv_model.hpp"
#include "hopf_observer.hpp"
#include "wave_disturbance.hpp"
#include "mpc_optim.hpp"
#include "actuator_saturation.hpp"
#include "telem_link.hpp"

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <fcntl.h>
#include <px4_platform_common/defines.h>
#include "trajectory.hpp"
/**
 * @file auv_hopf_main.cpp
 *
 * @brief AUV Hopf + Tube-MPC 主模块
 *
 * 当前版本定位：离线参考轨迹 + central/ancillary 两层 Tube-MPC 实现

 实物实验版本
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
    //TelemLink _telem{"/dev/ttyUSB0"};
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

    // PlotJuggler UDP 实时推流
    int _udp_sock{-1};
    struct sockaddr_in _pj_addr{};

    // mission-level 当前参考点（来自离线参考轨迹表，不直接等于第一层）
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
    float _sim_time{0.0f};

    // 离线参考轨迹索引
    int _ref_index{0};

    // 时域缓存
    float _mission_hor[auv::NX * auv::NHOR_MAX]{};
    float _mission_u_hor[auv::NU * auv::NHOR_MAX]{};
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
    float _preview_span_norm{0.0f};
    float _tube_tracking_error_norm{0.0f};
    float _tube_value_surrogate{0.0f};
    float _tube_value_delta{0.0f};
    float _central_state_violation_norm{0.0f};
    float _anc_state_violation_norm{0.0f};
    float _central_input_violation_norm{0.0f};
    float _anc_input_violation_norm{0.0f};
    float _physical_clip_norm{0.0f};
    bool _central_export_ok{false};
    bool _hopf_preview_ok{false};
    bool _central_bridge_ok{true};
    bool _anc_bridge_ok{true};

    // 调试输出
    uORB::Publication<auv_hopf_status_s> _status_pub{ORB_ID(auv_hopf_status)};

    // 内部帮助函数
    void reset_states();
    void build_mission_reference_from_table();
    void build_mission_reference_horizon_from_table();
    void build_disturbance_preview_horizon();
    void build_zero_disturbance_horizon();
    void build_central_rollout_horizon();
    void sync_nominal_from_plan_head();

    static void zero_horizon_buffer(float *buf, int total_len);
    static void shift_horizon_left(float *buf, int nhor, int dim);
    static void hold_last_stage(float *buf, int nhor, int dim);

    static void configure_mpc(MpcOptim &mpc,
                              const auv::AuvPlantParams &plant,
                              const auv::MpcWeightParams &weights,
                              const auv::MpcBoundsParams &bounds,
                              const auv::MpcStateBoundsParams &state_bounds,
                              const auv::MpcConfig &config);

    void clamp_state_horizon_to_bounds(float *x_hor,
                                       int nhor,
                                       const auv::MpcStateBoundsParams &bounds);
    void publish_debug(float dt);

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

    // [新增] 释放 UDP Socket 资源
    if (_udp_sock >= 0) {
        close(_udp_sock);
        _udp_sock = -1;
    }
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

float AuvHopf::wrap_pi(float angle)
{
    float res = fmodf(angle + M_PI_F, 2.0f * M_PI_F);
    if (res < 0.0f) {
        res += 2.0f * M_PI_F;
    }
    return res - M_PI_F;
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

void AuvHopf::zero_horizon_buffer(float *buf, int total_len)
{
    if (buf == nullptr || total_len <= 0) {
        return;
    }

    for (int i = 0; i < total_len; ++i) {
        buf[i] = 0.0f;
    }
}

void AuvHopf::shift_horizon_left(float *buf, int nhor, int dim)
{
    if (buf == nullptr || nhor < 2 || dim < 1) {
        return;
    }

    for (int k = 0; k < nhor - 1; ++k) {
        for (int i = 0; i < dim; ++i) {
            buf[k * dim + i] = buf[(k + 1) * dim + i];
        }
    }
}

void AuvHopf::hold_last_stage(float *buf, int nhor, int dim)
{
    if (buf == nullptr || nhor < 2 || dim < 1) {
        return;
    }

    for (int i = 0; i < dim; ++i) {
        buf[(nhor - 1) * dim + i] = buf[(nhor - 2) * dim + i];
    }
}

void AuvHopf::configure_mpc(MpcOptim &mpc,
                            const auv::AuvPlantParams &plant,
                            const auv::MpcWeightParams &weights,
                            const auv::MpcBoundsParams &bounds,
                            const auv::MpcStateBoundsParams &state_bounds,
                            const auv::MpcConfig &config)
{
    mpc.set_plant_params(plant);
    mpc.set_weight_params(weights);
    mpc.set_bounds_params(bounds);
    mpc.set_state_bounds_params(state_bounds);
    mpc.set_mpc_config(config);
    mpc.initialize(config.dt, config.Thor);
}

void AuvHopf::sync_nominal_from_plan_head()
{
    const int ksync = (_xref_nhor > 1) ? 1 : 0;

    for (int i = 0; i < auv::NX; ++i) {
        _z(i) = _zref_hor[ksync * auv::NX + i];
    }

    for (int i = 0; i < auv::NU; ++i) {
        _v_central(i) = _vref_hor[ksync * auv::NU + i];
    }
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
    _xref_nhor = 0;


    _disturbance_residual_norm = 0.0f;
    _preview_update_norm = 0.0f;
    _preview_span_norm = 0.0f;
    _tube_tracking_error_norm = 0.0f;
    _tube_value_surrogate = 0.0f;
    _tube_value_delta = 0.0f;
    _central_state_violation_norm = 0.0f;
    _anc_state_violation_norm = 0.0f;
    _central_input_violation_norm = 0.0f;
    _anc_input_violation_norm = 0.0f;
    _physical_clip_norm = 0.0f;
    _central_export_ok = false;
    _hopf_preview_ok = false;
    _central_bridge_ok = true;
    _anc_bridge_ok = true;

    zero_horizon_buffer(_mission_hor, auv::NX * auv::NHOR_MAX);
    zero_horizon_buffer(_mission_u_hor, auv::NU * auv::NHOR_MAX);
    zero_horizon_buffer(_zref_hor, auv::NX * auv::NHOR_MAX);
    zero_horizon_buffer(_vref_hor, auv::NU * auv::NHOR_MAX);
    zero_horizon_buffer(_dhat_hor, auv::NP * auv::NHOR_MAX);
    zero_horizon_buffer(_dhat_prev_hor, auv::NP * auv::NHOR_MAX);
    zero_horizon_buffer(_dzero_hor, auv::NP * auv::NHOR_MAX);

    _ref_index = _mission_ref_params.ref_start_index;
    if (_ref_index < 0) {
        _ref_index = 0;
    }
    if (_ref_index >= auv::TRAJ_LENGTH) {
        _ref_index = auv::TRAJ_LENGTH - 1;
    }

    // 初始状态对齐到离线参考轨迹起点
    _x(0) = auv::REF_TRAJECTORY[_ref_index][1]+2;
    _x(1) = auv::REF_TRAJECTORY[_ref_index][2]+2;
    _x(2) = auv::REF_TRAJECTORY[_ref_index][3];
    _x(3) = auv::REF_TRAJECTORY[_ref_index][4];
    _x(4) = auv::REF_TRAJECTORY[_ref_index][5];
    _x(5) = auv::REF_TRAJECTORY[_ref_index][6];

    _xref = _x;
    _z = _x;

    _wave.reset();
    _hopf.reset();
}

/**
 * @brief 从离线参考轨迹表读取当前 mission reference 点
 */
void AuvHopf::build_mission_reference_from_table()
{
    int idx = _ref_index;

    if (idx < 0) {
        idx = 0;
    }

    if (idx >= auv::TRAJ_LENGTH) {
        idx = auv::TRAJ_LENGTH - 1;
    }

    _xref(0) = auv::REF_TRAJECTORY[idx][1];
    _xref(1) = auv::REF_TRAJECTORY[idx][2];
    _xref(2) = auv::REF_TRAJECTORY[idx][3];
    _xref(3) = auv::REF_TRAJECTORY[idx][4];
    _xref(4) = auv::REF_TRAJECTORY[idx][5];
    _xref(5) = auv::REF_TRAJECTORY[idx][6];
}

void AuvHopf::build_mission_reference_horizon_from_table()
{
    _xref_nhor = (_central_mpc.horizon_length() > 0) ? _central_mpc.horizon_length() : 20;

    zero_horizon_buffer(_mission_hor, auv::NX * auv::NHOR_MAX);
    zero_horizon_buffer(_mission_u_hor, auv::NU * auv::NHOR_MAX);

    for (int k = 0; k < _xref_nhor && k < auv::NHOR_MAX; ++k) {
        int idx = _ref_index + k;

        if (idx >= auv::TRAJ_LENGTH) {
            idx = auv::TRAJ_LENGTH - 1; // 尾部保持
        }

        // REF_TRAJECTORY[idx] = {t, X, Y, psi, u, v, r, tau_u, tau_v, tau_r}
        for (int i = 0; i < auv::NX; ++i) {
            _mission_hor[k * auv::NX + i] = auv::REF_TRAJECTORY[idx][1 + i];
        }

        for (int i = 0; i < auv::NU; ++i) {
            _mission_u_hor[k * auv::NU + i] = auv::REF_TRAJECTORY[idx][7 + i];
        }
    }

    clamp_state_horizon_to_bounds(_mission_hor, _xref_nhor, _central_state_bounds);
}

void AuvHopf::build_disturbance_preview_horizon()
{
    const int nh = (_anc_mpc.horizon_length() > 1) ? _anc_mpc.horizon_length() : 2;
    const float dt_preview = (_anc_mpc.sample_time() > 0.0f)
        ? _anc_mpc.sample_time()
        : sanitize_dt(_central_mpc.sample_time());

    for (int k = 0; k < auv::NHOR_MAX; ++k) {
        for (int i = 0; i < auv::NP; ++i) {
            _dhat_hor[k * auv::NP + i] = 0.0f;
        }
    }

    // 优先使用 Hopf 内模的短时信号外推：
    //   d_hat[k+i|k] 不是“已知未来海况后的真扰动预测”，
    //   而是把当前已学得的振子内模状态复制出来，在短预测域上做自主前滚，
    //   从而得到供第二层使用的外生扰动预览序列。
    //   真实海浪—AUV 耦合中未被该内模捕获的部分，仍由 residual / tube 机制吸收。
    _hopf_preview_ok = _hopf.predict_horizon(_dhat_hor, nh, dt_preview);

    if (!_hopf_preview_ok) {
        for (int k = 0; k < nh && k < auv::NHOR_MAX; ++k) {
            _dhat_hor[k * auv::NP + 0] = _d_hat(0);
            _dhat_hor[k * auv::NP + 1] = _d_hat(1);
            _dhat_hor[k * auv::NP + 2] = _d_hat(2);
        }
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

    const float dspan_u = _dhat_hor[(nh - 1) * auv::NP + 0] - _dhat_hor[0];
    const float dspan_v = _dhat_hor[(nh - 1) * auv::NP + 1] - _dhat_hor[1];
    const float dspan_r = _dhat_hor[(nh - 1) * auv::NP + 2] - _dhat_hor[2];
    _preview_span_norm = sqrtf(dspan_u * dspan_u + dspan_v * dspan_v + dspan_r * dspan_r);

    for (int i = 0; i < auv::NP * auv::NHOR_MAX; ++i) {
        _dhat_prev_hor[i] = _dhat_hor[i];
    }
}

void AuvHopf::build_zero_disturbance_horizon()
{
    zero_horizon_buffer(_dzero_hor, auv::NP * auv::NHOR_MAX);
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

    zero_horizon_buffer(_zref_hor, auv::NX * auv::NHOR_MAX);
    zero_horizon_buffer(_vref_hor, auv::NU * auv::NHOR_MAX);

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

    configure_mpc(_central_mpc,
                  bundle.plant,
                  bundle.central_weights,
                  bundle.central_bounds,
                  bundle.central_state_bounds,
                  bundle.central_config);

    configure_mpc(_anc_mpc,
                  bundle.plant,
                  bundle.ancillary_weights,
                  bundle.ancillary_bounds,
                  bundle.ancillary_state_bounds,
                  bundle.ancillary_config);

    reset_states();


    _telem.init();

    PX4_INFO("auv_hopf started, profile=%s", _parameter_profile_name);

    // 初始化向 PlotJuggler 推流的 UDP Socket
    _udp_sock = socket(AF_INET, SOCK_DGRAM, 0);

    if (_udp_sock >= 0) {
        // [大师建议：非阻塞模式]
        // 即使网络出现丢包或缓冲区满，sendto 也会立即返回而不会卡住你的 AUV 控制主循环
        int flags = fcntl(_udp_sock, F_GETFL, 0);
        if (flags != -1) {
            fcntl(_udp_sock, F_SETFL, flags | O_NONBLOCK);
        }

        // 2. 配置目标地址（建议将端口和IP定义为常量或变量）
        _pj_addr.sin_family = AF_INET;
        _pj_addr.sin_port = htons(9870); // PlotJuggler 默认端口

        if (inet_pton(AF_INET, "127.0.0.1", &_pj_addr.sin_addr) <= 0) {
            PX4_ERR("Invalid PlotJuggler IP address");
        }

        PX4_INFO("PlotJuggler telemetry initialized on 127.0.0.1:9870");

        } else {
            // 如果创建失败，记录错误但不至于让整个模块崩溃（只是没法看图表）
            PX4_ERR("Failed to create PlotJuggler UDP socket: %d", errno);
        }


    return true;
}

void AuvHopf::run()
{
    if (!init()) {
        PX4_ERR("init failed");
        return;
    }

    float dt = (_anc_mpc.sample_time() > 0.0f) ? _anc_mpc.sample_time() : 0.02f;
    dt = sanitize_dt(dt);
    const hrt_abstime dt_us = static_cast<hrt_abstime>(dt * 1e6f);

    hrt_abstime next_run_time = hrt_absolute_time();

    while (!should_exit()) {

        next_run_time += dt_us;
        perf_begin(_perf_total_run);
        _sim_time += dt;

        if (_sim_time >= _experiment_cfg.sim_time_limit) {
            PX4_INFO("sim time limit reached: %.2f s", (double)_experiment_cfg.sim_time_limit);
            perf_end(_perf_total_run);
            break;
        }

        // 1) 当前 mission 参考点与第一层参考时域
        build_mission_reference_from_table();
        build_mission_reference_horizon_from_table();

        // 2) 真实扰动 / Hopf 学习 / 预览
        _d_true = _wave.step(_x, _sim_time);
        _d_hat  = _hopf.update(_x, _tau_prev, dt);
        _disturbance_residual_norm = input_norm(_d_true - _d_hat);

        build_disturbance_preview_horizon();
        build_zero_disturbance_horizon();

        // 3) 第一层：从当前 nominal state 出发，每拍求解 central MPC
        _central_mpc.set_reference_horizon_control_preview(
            _z,
            _mission_hor,
            _mission_u_hor,
            _dzero_hor,
            _xref_nhor);

        perf_begin(_perf_central_solve);
        _v_central = _central_mpc.solve();
        perf_end(_perf_central_solve);

        _central_bridge_ok = _central_mpc.debug().bridge_owner_match;
        _central_input_violation_norm =
            input_bound_violation_norm(_v_central, _central_mpc.bounds_params());

        build_central_rollout_horizon();
        sync_nominal_from_plan_head();

        // 4) 第二层：每拍跟踪 central path，并注入 Hopf 预览
        const int nh_ref = (_anc_mpc.horizon_length() < _xref_nhor)
                        ? _anc_mpc.horizon_length()
                        : _xref_nhor;

        _anc_mpc.set_reference_horizon_control_preview(
            _x,
            _zref_hor,
            _vref_hor,
            _dzero_hor ,
            nh_ref);

        perf_begin(_perf_anc_solve);
        _tau_anc = _anc_mpc.solve();
        perf_end(_perf_anc_solve);

        _anc_bridge_ok = _anc_mpc.debug().bridge_owner_match;
        _anc_input_violation_norm =
            input_bound_violation_norm(_tau_anc, _anc_mpc.bounds_params());

        // 5) 执行层限幅
        _tau_applied = _actuator_sat.apply(_tau_anc);
        _physical_clip_norm = input_norm(_actuator_sat.debug().tau_clip);


        char tx_buf[16];
        double stt = _tau_applied(2)*10.0f;
        snprintf(tx_buf, sizeof(tx_buf), "T%.1f\r\n", stt);
        _telem.send_string(tx_buf);

        _telem.update();
        float current_rpm = _telem.get_latest_rpm();

        // --- 安全检查：判断数据是否过期 ---
        uint64_t now_1 = hrt_absolute_time();
        // 如果超过 500000 微秒（0.5秒）没收到新数据，认为通信断开
        if (now_1 - _telem.get_last_rx_timestamp() > 500000) {
            PX4_WARN("MCU Telementry timeout! RPM data is stale.");
            // 这里可以触发 AUV 的安全保护模式，比如让电机停转，或者切换到纯位置估算模式
            current_rpm = 0.0f;
        }
        // PX4_INFO("velocity: %.2f m/s", (double)current_rpm);
        _tau_applied(2) = current_rpm / 10.0f; // 反向映射回控制输入
        //6) 推进真实系统
        _x = _model.step(_x, _tau_applied, _d_true, dt);
        _tau_prev = _tau_applied;

        // 7) 指标
        _tube_tracking_error_norm = tube_error_norm(_x, _z);
        _tube_value_delta =
            quadratic_tracking_value(_x, _z, _anc_mpc.weight_params()) - _tube_value_surrogate;
        _tube_value_surrogate =
            quadratic_tracking_value(_x, _z, _anc_mpc.weight_params());

        _central_state_violation_norm = state_bound_violation_norm(_z, _central_state_bounds);
        _anc_state_violation_norm = state_bound_violation_norm(_x, _anc_state_bounds);

        publish_debug(dt);

        // 8) 参考索引推进
        _ref_index += _mission_ref_params.ref_stride;
        if (_ref_index >= auv::TRAJ_LENGTH) {
            _ref_index = auv::TRAJ_LENGTH - 1;
        }
        perf_end(_perf_total_run);

        // 等待下一个采样拍
        const hrt_abstime now = hrt_absolute_time();
        if (next_run_time > now) {
            px4_usleep(next_run_time - now);
        } else {
            // 落后时直接追赶，不 sleep
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

Current version runs a cleaned two-layer architecture internally:
wave disturbance -> Hopf learner -> central-path MPC -> Hopf-augmented ancillary MPC -> physical saturation -> plant update

Key property of this version:
- mission reference is no longer treated as the central path itself
- central state z is generated by a dedicated nominal OCP
- ancillary layer explicitly tracks z* / v* and uses disturbance preview
- Hopf preview is generated by autonomous oscillator rollout instead of ZOH replication
- actuator_saturation is now the dedicated physical layer object
- runtime parameters are loaded through a named profile bundle: engineering_default
- debug stream is minimized to trajectory / disturbance-learning essentials
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

    // 1) 轨迹 / nominal / reference / 跟踪残差
    for (int i = 0; i < auv::NX; ++i) {
        msg.x[i]         = _x(i);
        msg.z[i]         = _z(i);
        msg.xref[i]      = _xref(i);
        msg.error_x_z[i] = _x(i) - _z(i);
    }

    // 2) 干扰：真实 / 学习，以及实际控制量
    for (int i = 0; i < auv::NU; ++i) {
        msg.d_true[i] = _d_true(i);
        msg.d_hat[i]  = _d_hat(i);
        msg.d_res[i]  = _d_true(i) - _d_hat(i);

        // 新增：把实际控制量存入 uORB 消息，以便记录到 .ulg 日志中
        msg.tau_applied[i] = _tau_applied(i);
    }
    // 3) 仅保留跟踪残差范数（uORB 侧最小化）
    msg.tube_tracking_error_norm = _tube_tracking_error_norm;

    _status_pub.publish(msg);

    // 4) PlotJuggler：轨迹 + 干扰 + 每通道残差
    if (_udp_sock >= 0) {
        char json_buf[1024];

        const float ru = _d_true(0) - _d_hat(0);
        const float rv = _d_true(1) - _d_hat(1);
        const float rr = _d_true(2) - _d_hat(2);

        const int len = snprintf(
            json_buf, sizeof(json_buf),
            "{"
            "\"traj\":{"
                "\"x_actual\":%.3f, \"y_actual\":%.3f, "
                "\"x_nominal\":%.3f, \"y_nominal\":%.3f, "
                "\"x_ref\":%.3f, \"y_ref\":%.3f"
            "}, "
            "\"track_res\":{"
                "\"ex\":%.3f, \"ey\":%.3f, \"norm\":%.3f"
            "}, "
            "\"dist_surge\":{"
                "\"true\":%.3f, \"hat\":%.3f, \"res\":%.3f"
            "}, "
            "\"dist_sway\":{"
                "\"true\":%.3f, \"hat\":%.3f, \"res\":%.3f"
            "}, "
            "\"dist_yaw\":{"
                "\"true\":%.3f, \"hat\":%.3f, \"res\":%.3f"
            "}, "
            "\"tau_applied\":{"
                "\"u\":%.3f, \"v\":%.3f, \"r\":%.3f"
            "}"
            "}",
            (double)_x(0),    (double)_x(1),
            (double)_z(0),    (double)_z(1),
            (double)_xref(0), (double)_xref(1),

            (double)(_x(0) - _z(0)),
            (double)(_x(1) - _z(1)),
            (double)_tube_tracking_error_norm,

            (double)_d_true(0), (double)_d_hat(0), (double)ru,
            (double)_d_true(1), (double)_d_hat(1), (double)rv,
            (double)_d_true(2), (double)_d_hat(2), (double)rr,

            (double)_tau_applied(0), (double)_tau_applied(1), (double)_tau_applied(2)
        );

        if (len > 0 && len < (int)sizeof(json_buf)) {
            sendto(_udp_sock, json_buf, len, 0,
                   (struct sockaddr *)&_pj_addr, sizeof(_pj_addr));
        } else {
            PX4_WARN("JSON buffer truncated or error!");
        }
    }
}

ModuleBase::Descriptor AuvHopf::desc{
    AuvHopf::task_spawn, AuvHopf::custom_command, AuvHopf::print_usage
};

extern "C" __EXPORT int auv_hopf_main(int argc, char *argv[])
{
    return ModuleBase::main(AuvHopf::desc, argc, argv);
}
