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
#include "thruster_allocation.hpp"
#include "telem_link.hpp"

/**
 * @file auv_hopf_main.cpp
 *
 * @brief AUV Hopf + Tube-MPC 主模块
 *
 * 当前版本定位：
 * 1. 先把学术代码主链在 PX4 模块内部完整跑通
 * 2. 使用内部状态 x / z 做闭环仿真
 * 3. 后续再根据需要接入真实 PX4 topic
 *
 * 主要数据流：
 *   d_true      = wave_disturbance.step(x, t)
 *   d_hat       = hopf_observer.update(x, tau_prev, dt)
 *   tau_nom     = mpc_optim.solve(z, xref_hor, d_hat)
 *   tau_applied = thruster_allocation.compute(x, z, d_hat, tau_nom)
 *   x_next      = model.step(x, tau_applied, d_true, dt)
 *   z_next      = model.step(z, tau_nom, d_hat, dt)
 *
 * 说明：
 * - x 为真实状态（actual state）
 * - z 为名义状态（nominal state）
 * - 二者在 reset 时可相同，但运行中必须独立推进
 * - 这样 ancillary / tube 误差 e = x - z 才有真实控制意义
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
    MpcOptim _mpc;
    ThrusterAllocation _alloc;
    TelemLink _telem{"/dev/ttyS1"}; // <--- 实例化串口模块

    //----------------------------------------------
    // 性能计数器
    //----------------------------------------------
    perf_counter_t _perf_mpc_solve{perf_alloc(PC_ELAPSED, "auv_hopf_mpc")};
    perf_counter_t _perf_total_run{perf_alloc(PC_ELAPSED, "auv_hopf_total")};

    // --------------------------------------------------
    // 状态变量
    // x : 真实状态 actual state
    // z : 名义状态 nominal state
    // --------------------------------------------------
    auv::StateVec _x{};
    auv::StateVec _z{};

    // 当前参考点
    auv::StateVec _xref{};

    // 当前真实扰动 / 学习估计 / 控制量
    auv::InputVec _d_true{};
    auv::InputVec _d_hat{};
    auv::InputVec _tau_nom{};
    auv::InputVec _tau_applied{};
    auv::InputVec _tau_prev{};

    // 调度时间
    hrt_abstime _last_run{0};
    float _sim_time{0.0f};

    // MPC 参考时域缓存
    float _xref_hor[auv::NX * auv::NHOR_MAX]{};
    int _xref_nhor{0};

    // 调试输出
    uORB::Publication<auv_hopf_status_s> _status_pub{ORB_ID(auv_hopf_status)};

    // 内部帮助函数
    void reset_states();
    void build_reference(float t);
    void build_reference_horizon();
    void publish_debug(float dt);

    static float sanitize_dt(float dt);
};

// 构造函数
AuvHopf::AuvHopf() :
    ModuleParams(nullptr)
{
    reset_states();
}
// 析构函数
AuvHopf::~AuvHopf()
{
    perf_free(_perf_mpc_solve);
    perf_free(_perf_total_run);
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
    _tau_nom.zero();
    _tau_applied.zero();
    _tau_prev.zero();

    _sim_time = 0.0f;
    _last_run = 0;
    _xref_nhor = 0;

    // -----------------------------
    // 初始真实状态
    // -----------------------------
    _x(0) = 0.0f;   // X
    _x(1) = 0.0f;   // Y
    _x(2) = 0.0f;   // psi
    _x(3) = 0.0f;   // u
    _x(4) = 0.0f;   // v
    _x(5) = 0.0f;   // r

    // -----------------------------
    // 初始名义状态
    //
    // reset 时允许 z(0)=x(0)
    // 但运行后 z 必须独立递推，不能再被 x 覆盖
    // -----------------------------
    _z = _x;

    for (int i = 0; i < auv::NX * auv::NHOR_MAX; ++i) {
        _xref_hor[i] = 0.0f;
    }

    _wave.reset();
    _hopf.reset();
    _mpc.reset();
}

/**
 * @brief 构造当前参考轨迹点
 *
 * 示例使用平面圆轨迹：
 *   Xref = R cos(wt)
 *   Yref = R sin(wt)
 *   psiref = atan2(Ydot, Xdot)
 */
void AuvHopf::build_reference(float t)
{
    const float R = 30.0f;
    const float w = 0.02f;

    const float xr = R * cosf(w * t);
    const float yr = R * sinf(w * t);

    const float xdot_r = -R * w * sinf(w * t);
    const float ydot_r =  R * w * cosf(w * t);

    const float psi_r = atan2f(ydot_r, xdot_r);
    const float u_r = R*w;
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
 * @brief 构造 MPC 参考时域
 *
 * 按平铺方式组织：
 *   [x0(0..5), x1(0..5), ..., xN(0..5)]
 */
void AuvHopf::build_reference_horizon()
{
    _xref_nhor = (_mpc.horizon_length() > 0) ? _mpc.horizon_length() : 20;
    const float dt = (_mpc.sample_time() > 0.0f) ? _mpc.sample_time() : 0.02f;

    for (int k = 0; k < auv::NHOR_MAX; ++k) {
        for (int i = 0; i < auv::NX; ++i) {
            _xref_hor[k * auv::NX + i] = 0.0f;
        }
    }

    for (int k = 0; k < _xref_nhor && k < auv::NHOR_MAX; ++k) {
        const float tk = _sim_time + static_cast<float>(k) * dt;

        const float R = 30.0f;
        const float w = 0.02f;

        const float xr = R * cosf(w * tk);
        const float yr = R * sinf(w * tk);

        const float xdot_r = -R * w * sinf(w * tk);
        const float ydot_r =  R * w * cosf(w * tk);

        const float psi_r = atan2f(ydot_r, xdot_r);
        const float u_r = R * w;
        const float v_r = 0.0f;
        const float r_r = w;

        _xref_hor[k * auv::NX + 0] = xr;
        _xref_hor[k * auv::NX + 1] = yr;
        _xref_hor[k * auv::NX + 2] = psi_r;
        _xref_hor[k * auv::NX + 3] = u_r;
        _xref_hor[k * auv::NX + 4] = v_r;
        _xref_hor[k * auv::NX + 5] = r_r;
    }
}

bool AuvHopf::init()
{
    const auto plant = auv::default_plant_params();
    const auto hopf_cfg = auv::default_hopf_config();
    const auto wave_env = auv::default_wave_env_params();
    const auto wave_model = auv::default_wave_force_model_params();
    const auto tube = auv::default_tube_feedback_params();
    const auto mpc_w = auv::default_mpc_weight_params();
    const auto mpc_b = auv::default_mpc_bounds_params();
    const auto mpc_c = auv::default_mpc_config();

    _model.set_plant_params(plant);

    _hopf.set_plant_params(plant);
    _hopf.set_hopf_config(hopf_cfg);

    _wave.set_wave_env(wave_env);
    _wave.set_model_params(wave_model);

    _alloc.set_plant_params(plant);
    _alloc.set_params(tube);

    _mpc.set_plant_params(plant);
    _mpc.set_weight_params(mpc_w);
    _mpc.set_bounds_params(mpc_b);
    _mpc.set_mpc_config(mpc_c);
    _mpc.initialize(mpc_c.dt, mpc_c.Thor);

    reset_states();

    _telem.init();

    PX4_INFO("auv_hopf started");
    return true;
}

void AuvHopf::run()
{
    if (!init()) {
        PX4_ERR("init failed");
        return;
    }

    float dt = (_mpc.sample_time() > 0.0f) ? _mpc.sample_time() : 0.02f;
    dt = sanitize_dt(dt);
    const hrt_abstime dt_us = static_cast<hrt_abstime>(dt * 1e6f);

    //锚定系统绝对时间作为循环起点，避免初始 dt 过大
    hrt_abstime next_run_time = hrt_absolute_time();
    // TelemData current_telem_data{};

    while (!should_exit()) {

        next_run_time += dt_us;
        //---------------------------------
        // --- 算法执行块开始 性能计数器：总运行时间---
        perf_begin(_perf_total_run);
        _sim_time += dt;

        // 尝试提取最新数据
        // if (_telem.get_latest_data(current_telem_data)) {
        //     // 如果拿到了新数据，就在这里把数据喂给你的算法！
        //     // 比如： _x(3) = current_telem_data.rpm * 0.1f;
        //     PX4_INFO("Got fresh data! RPM = %.2f", (double)current_telem_data.rpm);
        // } else {
        //     // 如果没拿到新数据（串口没发或者频率没对上）
        //     // 算法会静默使用上一拍的老数据继续算，保证绝对安全！
        // }

        // --------------------------------------------------
        // 1) 构造当前参考与参考时域
        build_reference(_sim_time);
        build_reference_horizon();

        // --------------------------------------------------
        // 2) 真实外扰作用
        _d_true = _wave.step(_x, _sim_time);

        // --------------------------------------------------
        // 3) Hopf 在线扰动学习
        _d_hat = _hopf.update(_x, _tau_prev, dt);

        // --------------------------------------------------
        // 4) 用 nominal state z 设置 MPC 初值和时域参考
        auv::ParamVec d_zero{};
        d_zero.zero();
        _mpc.set_reference_horizon(_z, _xref_hor, _xref_nhor, d_zero);

        // --------------------------------------------------
        // 5) 求 nominal control
        // 性能计数器：MPC 求解时间
        perf_begin(_perf_mpc_solve); // 核心：单独记录 MPC 求解耗时
        _tau_nom = _mpc.solve();
        perf_end(_perf_mpc_solve);

        // --------------------------------------------------
        // 6) ancillary / allocation
        _tau_applied = _alloc.compute(_x, _z, _d_hat, _tau_nom);

        // --------------------------------------------------
        // 7) 真实状态推进：actual state x
        _x = _model.step(_x, _tau_applied, _d_true, dt);

        // --------------------------------------------------
        // 8) 名义状态推进：nominal state z
        _z = _model.step(_z, _tau_nom, d_zero, dt);

        // --------------------------------------------------
        // 9) 更新上一拍实际执行控制
        _tau_prev = _tau_applied;

        // 如果算出来是 6.82，就会拼成 "T6.8\r\n" 发送出去
        // char tx_buf[16];
        // //snprintf(tx_buf, sizeof(tx_buf), "T%.1f\r\n", (double)_tau_applied(0));
        // double stt = 6.8;
        // snprintf(tx_buf, sizeof(tx_buf), "T%.1f\r\n", stt);
        // _telem.send_string(tx_buf);
        // --------------------------------------------------
        // 10) 发布调试信息
        publish_debug(dt);

        // --------------------------------------------------
        // 算法执行结束，性能计数器：总运行时间
        perf_end(_perf_total_run);
        // --------------------------------------------------
        // 4. 时序控制：计算剩余裕量并让出 CPU
        hrt_abstime now = hrt_absolute_time();

        if (_sim_time >= 300.0f) {
                    PX4_INFO("时间到！已运行 300 秒，算法自动安全退出。");
                    request_stop(); // 优雅地触发 should_exit() = true
                }

        if (now < next_run_time) {
            px4_usleep(next_run_time - now);
        } else {
            PX4_WARN("Deadline Missed! dt overrun.");
            next_run_time = now; // 重置时间锚点
        }
    }
}

int AuvHopf::print_status()
{
    PX4_INFO("running");
    PX4_INFO("t = %.3f", (double)_sim_time);

    PX4_INFO("x = [%.3f %.3f %.3f %.3f %.3f %.3f]",
             (double)_x(0), (double)_x(1), (double)_x(2),
             (double)_x(3), (double)_x(4), (double)_x(5));

    PX4_INFO("z = [%.3f %.3f %.3f %.3f %.3f %.3f]",
             (double)_z(0), (double)_z(1), (double)_z(2),
             (double)_z(3), (double)_z(4), (double)_z(5));

    PX4_INFO("xref = [%.3f %.3f %.3f %.3f %.3f %.3f]",
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

    PX4_INFO("tau_nom = [%.3f %.3f %.3f]",
             (double)_tau_nom(0), (double)_tau_nom(1), (double)_tau_nom(2));

    PX4_INFO("tau_applied = [%.3f %.3f %.3f]",
             (double)_tau_applied(0), (double)_tau_applied(1), (double)_tau_applied(2));

    PX4_INFO("xref_nhor = %d", _xref_nhor);

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

    // 将实例注册到系统描述符中
    desc.object.store(instance);

    // 开始执行核心硬实时大循环 (run 内部会自动调用 init)
    instance->run();

    // 只有当收到 stop 命令，退出 while(!should_exit()) 后，才会执行到这里
    delete instance;
    desc.object.store(nullptr);
    desc.task_id = -1;
    return 0;
}

// 替换原来的 task_spawn
int AuvHopf::task_spawn(int argc, char *argv[])
{
    // 向系统申请独立线程，并将程序入口指向 trampoline
    int task_id = px4_task_spawn_cmd("auv_hopf",
                                     SCHED_DEFAULT,
                                     SCHED_PRIORITY_ATTITUDE_CONTROL,
                                     4096,
                                     (px4_main_t)&AuvHopf::trampoline, // 指向跳板函数
                                     (char *const *)argv);

    if (task_id < 0) {
        return -errno; // 启动失败
    }

    // 🚨 把系统分配给你的真实特权线程号存起来
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

Current version runs the academic closed-loop chain internally:
wave disturbance -> Hopf learner -> MPC -> ancillary allocation -> plant update

Key property of this version:
- actual state x and nominal state z are propagated independently
- z is no longer overwritten by x
- tube / ancillary feedback therefore has actual closed-loop meaning
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

    // 状态与误差拷贝
    for (int i = 0; i < auv::NX; ++i) {
        msg.x[i]         = _x(i);
        msg.z[i]         = _z(i);
        msg.xref[i]      = _xref(i);
        msg.error_x_z[i] = _x(i) - _z(i);
    }

    // 3维向量拷贝
    for (int i = 0; i < auv::NU; ++i) {
        msg.d_true[i]      = _d_true(i);
        msg.d_hat[i]       = _d_hat(i);
        msg.tau_nom[i]     = _tau_nom(i);
        msg.tau_applied[i] = _tau_applied(i);
    }

    // 内部调试信息
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

    const auto &alloc_dbg = _alloc.debug();
    for (int i = 0; i < 3; ++i) {
        msg.tau_fb[i]        = alloc_dbg.tau_fb(i);
        msg.tau_ancillary[i] = alloc_dbg.tau_ancillary(i);
        msg.tau_raw[i]       = alloc_dbg.tau_raw(i);
    }

    msg.xref_nhor = static_cast<uint8_t>(_xref_nhor);

    // 发布！
    _status_pub.publish(msg);
}


ModuleBase::Descriptor AuvHopf::desc{AuvHopf::task_spawn, AuvHopf::custom_command, AuvHopf::print_usage};

extern "C" __EXPORT int auv_hopf_main(int argc, char *argv[])
{
    // 正确的写法：调用 ModuleBase 的 main，并把我们自己的名片 (desc) 传进去！
    return ModuleBase::main(AuvHopf::desc, argc, argv);
}

