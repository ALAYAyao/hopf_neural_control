#pragma once

/**
 * @file auv_common.hpp
 *
 * @brief AUV Hopf + Tube-MPC 项目的公共类型、参数与默认配置
 *
 * 本文件是整个项目的统一公共入口，负责集中定义：
 * 1. 基本维度与公共向量类型
 * 2. AUV 名义动力学参数
 * 3. Hopf 学习器参数
 * 4. 波浪环境与波浪力模型参数
 * 5. Tube / ancillary 反馈参数
 * 6. MPC / GRAMPC 参数、约束与默认配置
 *
 * 状态:
 *   x = [X, Y, psi, u, v, r]
 *
 * 控制 / 扰动:
 *   tau = [tau_u, tau_v, tau_r]
 *   d   = [d_u, d_v, d_r]
 */

#include <matrix/matrix/math.hpp>
#include <stdint.h>

namespace auv
{

// ==================================================
// 基本维度定义
// ==================================================

static constexpr int NX = 6;      // [X, Y, psi, u, v, r]
static constexpr int NU = 3;      // [tau_u, tau_v, tau_r]
static constexpr int NP = 3;      // [d_u, d_v, d_r]
static constexpr int N_OSC = 4;   // 每个通道 Hopf 振子数

// 来源版 wave_disturbance 使用 4 个频率点、4 个方向点
static constexpr int N_WAVE_FREQ = 4;
static constexpr int N_WAVE_DIR  = 4;

// MPC 预测域参考最大长度
static constexpr int NHOR_MAX = 64;

using StateVec = matrix::Vector<float, NX>;
using InputVec = matrix::Vector<float, NU>;
using ParamVec = matrix::Vector<float, NP>;

// ==================================================
// AUV 名义动力学参数
// ==================================================

/**
 * @brief AUV 名义动力学参数
 *
 * 对应 3-DOF 平面模型中的质量/惯量与线性阻尼。
 */
struct AuvPlantParams {
    float m_u;
    float m_v;
    float Iz;

    float Du;
    float Dv;
    float Dr;
};

// ==================================================
// Hopf 学习器参数
// ==================================================

struct HopfParams {
    float gamma;
    float mu;
    float epsF;
    float epsW;
    float eta;
    float lambda_alpha;
};

struct HopfChannelInit {
    float xi0[N_OSC];
    float yi0[N_OSC];
    float wi0[N_OSC];
};

struct HopfConfig {
    HopfParams surge;
    HopfParams sway;
    HopfParams yaw;

    HopfChannelInit init_surge;
    HopfChannelInit init_sway;
    HopfChannelInit init_yaw;
};

// ==================================================
// 波浪环境参数
// ==================================================

/**
 * @brief 来源版 wave_env 参数
 *
 * 对应原始代码中的:
 *   wave_env.omegaMax
 *   wave_env.Omega
 *   wave_env.S
 *   wave_env.mu
 *   wave_env.d
 *   wave_env.Epsilon
 */
struct WaveEnvParams {
    float omegaMax;

    // 频率采样点
    float Omega[N_WAVE_FREQ];

    // 一维频谱值（与频率点对应）
    float S[N_WAVE_FREQ];

    // 方向采样点
    float mu[N_WAVE_DIR];

    // 水深
    float d;

    // 随机相位矩阵 Epsilon(n,m)
    float Epsilon[N_WAVE_FREQ][N_WAVE_DIR];
};

// ==================================================
// 波浪力模型参数
// ==================================================

/**
 * @brief waveforce_SIM 中使用的 AUV 物理与数值参数
 */
struct WaveForceModelParams {
    float rho;           // 海水密度
    float L;             // AUV 长度
    float D_auv;         // AUV 直径
    float cd;            // 阻力系数
    float cm1;           // surge / 轴向附加惯性相关系数
    float cm2;           // 横向附加惯性相关系数
    float g;             // 重力加速度
    float h;             // 沿艇体离散步长
    float scale_factor;  // wave_disturbance 外层输出缩放
};

// ==================================================
// 约束参数
// ==================================================

/**
 * @brief AUV 约束参数
 *
 * 当前主要作为统一的约束配置入口。
 * 后续若需要显式状态/输入约束建模，可继续扩展。
 */
struct AuvConstraintParams {
    // 输入约束
    float tau_u_min;
    float tau_u_max;

    float tau_v_min;
    float tau_v_max;

    float tau_r_min;
    float tau_r_max;

    // 状态约束
    float psi_min;
    float psi_max;

    float u_min;
    float u_max;

    float v_min;
    float v_max;

    float r_min;
    float r_max;
};

// ==================================================
// Tube / ancillary feedback 参数
// ==================================================

/**
 * @brief Tube-MPC 辅助反馈与执行层限幅参数
 *
 * 注意：
 * 1. 这里只保留 ancillary / allocation 层真正需要的反馈增益与物理限幅；
 * 2. 质量/惯量参数不再重复保存，统一由 AuvPlantParams 提供；
 * 3. 这样可以避免 M11/M22/M33 与 m_u/m_v/Iz 长期漂移。
 */
struct TubeFeedbackParams {
    // 位置 / 姿态相关比例反馈
    float Kp11;
    float Kp22;
    float Kp33;

    // 速度反馈
    float Kv11;
    float Kv22;
    float Kv33;

    // 执行层物理限幅
    float tau_u_max_phys;
    float tau_v_max_phys;
    float tau_r_max_phys;
};


/**
 * @brief 执行层物理限幅参数
 *
 * 与 TubeFeedbackParams 拆开保存，避免“反馈器参数”和“执行器硬限幅”继续绑在一起。
 */
struct PhysicalActuatorLimits {
    float tau_u_min;
    float tau_u_max;

    float tau_v_min;
    float tau_v_max;

    float tau_r_min;
    float tau_r_max;
};

/**
 * @brief 任务层参考轨迹参数
 */
struct MissionReferenceParams {
    float circle_radius;
    float circle_omega;
};

/**
 * @brief 运行实验配置
 */
struct ExperimentConfig {
    float sim_time_limit;
};

/**
 * @brief 参数口径分组
 *
 * EngineeringDefault: 面向当前工程原型迭代。
 * PaperReproduction: 面向论文复现实验口径。
 */
enum class ParameterProfile : uint8_t {
    EngineeringDefault = 0,
    PaperReproduction = 1,
};

/**
 * @brief 主模块一次性装载的参数包
 *
 * 这样可以显式区分“论文复现参数组”和“工程默认参数组”，
 * 避免后续主模块继续从不同默认工厂里零散拼装，导致口径漂移。
 */
struct RuntimeProfileBundle {
    const char *profile_name;

    AuvPlantParams plant;
    HopfConfig hopf;
    WaveEnvParams wave_env;
    WaveForceModelParams wave_model;

    TubeFeedbackParams tube_feedback;
    PhysicalActuatorLimits physical_limits;
    MissionReferenceParams mission_ref;
    ExperimentConfig experiment;

    MpcWeightParams central_weights;
    MpcBoundsParams central_bounds;
    MpcStateBoundsParams central_state_bounds;
    MpcConfig central_config;

    MpcWeightParams ancillary_weights;
    MpcBoundsParams ancillary_bounds;
    MpcStateBoundsParams ancillary_state_bounds;
    MpcConfig ancillary_config;
};

// ==================================================
// MPC / GRAMPC 参数
// ==================================================

/**
 * @brief MPC 代价权重
 *
 * 阶段代价:
 *   l = 0.5 * (dx'Qdx + du'Rdu)
 *
 * 终端代价:
 *   V = 0.5 * dx'Pdx
 */
struct MpcWeightParams {
    float Q[NX];
    float R[NU];
    float P[NX];
};

/**
 * @brief MPC 输入/参数边界
 *
 * 当前 p = [d_u, d_v, d_r] 作为 GRAMPC 的已知参数输入。
 */
struct MpcBoundsParams {
    float umin[NU];
    float umax[NU];

    float pmin[NP];
    float pmax[NP];
};

/**
 * @brief MPC 软状态边界与罚函数配置
 *
 * 说明：
 * 1. 当前位置版本先把状态边界显式送进优化代价，而不是只留在文档里；
 * 2. 由于当前 GRAMPC 接口仍采用无显式 h/g 约束配置，这里先使用 soft-constraint 方式；
 * 3. 这样可以先把 nominal / ancillary 的状态安全边界真正区分出来，
 *    后续若切到显式不等式约束，可在保持参数口径不变的前提下升级实现。
 */
struct MpcStateBoundsParams {
    bool enable_soft_bounds;

    float xmin[NX];
    float xmax[NX];

    // 对越界量 0.5 * w_i * viol_i^2 计入阶段/终端代价
    float penalty[NX];
};

/**
 * @brief MPC 求解器配置
 */
struct MpcConfig {
    float dt;
    float Thor;

    int Nhor;
    int MaxGradIter;
    int MaxMultIter;

    bool shift_control;
    bool use_terminal_cost;
    bool use_integral_cost;
};

// ==================================================
// 默认工厂函数
// ==================================================

inline AuvPlantParams default_plant_params()
{
    AuvPlantParams p{};

    p.m_u = 31.8381f;
    p.m_v = 60.4368f;
    p.Iz  = 7.4501f;

    p.Du  = 1.5919f;
    p.Dv  = 3.0218f;
    p.Dr  = 7.4501f;

    return p;
}

inline HopfConfig default_hopf_config()
{
    HopfConfig cfg{};

    cfg.surge.gamma = 1.0f;
    cfg.surge.mu = 1.0f;
    cfg.surge.epsF = 0.7f;
    cfg.surge.epsW = 0.5f;
    cfg.surge.eta = 1.0f;
    cfg.surge.lambda_alpha = 0.0f;

    cfg.sway.gamma = 1.0f;
    cfg.sway.mu = 1.0f;
    cfg.sway.epsF = 0.7f;
    cfg.sway.epsW = 0.55f;
    cfg.sway.eta = 0.5f;
    cfg.sway.lambda_alpha = 0.0f;

    cfg.yaw.gamma = 3.0f;
    cfg.yaw.mu = 1.0f;
    cfg.yaw.epsF = 0.3f;
    cfg.yaw.epsW = 0.3f;
    cfg.yaw.eta = 3.0f;
    cfg.yaw.lambda_alpha = 0.06f;

    const float w0[N_OSC]  = {0.5f, 1.3333f, 2.1667f, 3.0f};

    const float xu0[N_OSC] = { 0.0291f,  0.0083f,  0.0138f, -0.0106f};
    const float yu0[N_OSC] = {-0.0047f, -0.0027f,  0.0110f, -0.0028f};

    const float xv0[N_OSC] = { 0.0070f, -0.0205f, -0.0035f, -0.0082f};
    const float yv0[N_OSC] = {-0.0158f,  0.0051f,  0.0028f,  0.0003f};

    const float xr0[N_OSC] = {-0.0133f,  0.0113f,  0.0035f, -0.0030f};
    const float yr0[N_OSC] = { 0.0002f, -0.0026f, -0.0175f, -0.0029f};

    for (int i = 0; i < N_OSC; ++i) {
        cfg.init_surge.xi0[i] = xu0[i];
        cfg.init_surge.yi0[i] = yu0[i];
        cfg.init_surge.wi0[i] = w0[i];

        cfg.init_sway.xi0[i] = xv0[i];
        cfg.init_sway.yi0[i] = yv0[i];
        cfg.init_sway.wi0[i] = w0[i];

        cfg.init_yaw.xi0[i] = xr0[i];
        cfg.init_yaw.yi0[i] = yr0[i];
        cfg.init_yaw.wi0[i] = w0[i];
    }

    return cfg;
}

inline WaveEnvParams default_wave_env_params()
{
    WaveEnvParams w{};

    w.omegaMax = 3.0f;

    w.Omega[0] = 0.5000f;
    w.Omega[1] = 1.3333f;
    w.Omega[2] = 2.1667f;
    w.Omega[3] = 3.0000f;

    w.S[0] = 0.0f;
    w.S[1] = 0.0058f;
    w.S[2] = 0.0099f;
    w.S[3] = 0.0028f;

    w.mu[0] = -1.5708f;
    w.mu[1] = -0.5236f;
    w.mu[2] =  0.5236f;
    w.mu[3] =  1.5708f;

    w.d = 3.0f;

    const float eps[N_WAVE_FREQ][N_WAVE_DIR] = {
        {4.6321f, 1.9052f, 4.0146f, 2.6159f},
        {1.1387f, 4.9017f, 5.6118f, 4.6482f},
        {2.6782f, 3.3738f, 0.3813f, 5.6106f},
        {0.6168f, 4.8322f, 1.1044f, 0.1624f}
    };

    for (int i = 0; i < N_WAVE_FREQ; ++i) {
        for (int j = 0; j < N_WAVE_DIR; ++j) {
            w.Epsilon[i][j] = eps[i][j];
        }
    }

    return w;
}

inline WaveForceModelParams default_wave_force_model_params()
{
    WaveForceModelParams p{};

    p.rho = 1026.0f;
    p.L = 1.6f;
    p.D_auv = 0.19f;
    p.cd = 0.42f;
    p.cm1 = 1.03f;
    p.cm2 = 1.9652f;
    p.g = 9.7803267714f;
    p.h = 0.1f;
    p.scale_factor = 0.5f;

    return p;
}

inline AuvConstraintParams default_constraint_params()
{
    AuvConstraintParams c{};

    c.tau_u_min = -50.0f;
    c.tau_u_max =  50.0f;

    c.tau_v_min = -50.0f;
    c.tau_v_max =  50.0f;

    c.tau_r_min = -10.0f;
    c.tau_r_max =  10.0f;

    c.psi_min = -3.1415926f;
    c.psi_max =  3.1415926f;

    c.u_min = -3.0f;
    c.u_max =  3.0f;

    c.v_min = -3.0f;
    c.v_max =  3.0f;

    c.r_min = -2.0f;
    c.r_max =  2.0f;

    return c;
}

inline TubeFeedbackParams default_tube_feedback_params()
{
    TubeFeedbackParams p{};

    // ancillary_and_allocation 中的 Kp
    p.Kp11 = 8.0f;
    p.Kp22 = 8.0f;
    p.Kp33 = 4.0f;

    // ancillary_and_allocation 中的 Kv
    p.Kv11 = 15.0f;
    p.Kv22 = 15.0f;
    p.Kv33 = 8.0f;

    // 物理限幅
    p.tau_u_max_phys = 50.0f;
    p.tau_v_max_phys = 50.0f;
    p.tau_r_max_phys = 10.0f;

    return p;
}

inline PhysicalActuatorLimits default_physical_actuator_limits()
{
    PhysicalActuatorLimits p{};

    p.tau_u_min = -50.0f;
    p.tau_u_max =  50.0f;

    p.tau_v_min = -50.0f;
    p.tau_v_max =  50.0f;

    p.tau_r_min = -10.0f;
    p.tau_r_max =  10.0f;

    return p;
}

inline PhysicalActuatorLimits paper_reproduction_physical_actuator_limits()
{
    PhysicalActuatorLimits p{};

    p.tau_u_min = -50.0f;
    p.tau_u_max =  50.0f;

    p.tau_v_min = -50.0f;
    p.tau_v_max =  50.0f;

    // 论文当前给出的物理 yaw 力矩限幅为 ±20
    p.tau_r_min = -20.0f;
    p.tau_r_max =  20.0f;

    return p;
}

inline MissionReferenceParams default_mission_reference_params()
{
    MissionReferenceParams p{};
    p.circle_radius = 30.0f;
    p.circle_omega  = 0.02f;
    return p;
}

inline MissionReferenceParams paper_reproduction_mission_reference_params()
{
    MissionReferenceParams p{};
    p.circle_radius = 20.0f;
    p.circle_omega  = 0.05f;
    return p;
}

inline ExperimentConfig default_experiment_config()
{
    ExperimentConfig cfg{};
    cfg.sim_time_limit = 300.0f;
    return cfg;
}

inline ExperimentConfig paper_reproduction_experiment_config()
{
    ExperimentConfig cfg{};
    cfg.sim_time_limit = 60.0f;
    return cfg;
}

inline AuvPlantParams paper_reproduction_plant_params()
{
    AuvPlantParams p{};

    p.m_u = 70.0f;
    p.m_v = 75.0f;
    p.Iz  = 12.501f;

    p.Du  = 15.919f;
    p.Dv  = 20.218f;
    p.Dr  = 10.4501f;

    return p;
}

inline MpcWeightParams default_central_mpc_weight_params()
{
    MpcWeightParams w{};

    const float q[NX] = {2000.0f, 2000.0f, 1000.0f, 0.2f, 0.2f, 0.2f};
    const float r[NU] = {0.1f, 0.1f, 0.1f};
    const float p[NX] = {2000.0f, 2000.0f, 1000.0f, 0.2f, 0.2f, 0.2f};

    for (int i = 0; i < NX; ++i) {
        w.Q[i] = q[i];
        w.P[i] = p[i];
    }

    for (int i = 0; i < NU; ++i) {
        w.R[i] = r[i];
    }

    return w;
}

inline MpcWeightParams default_ancillary_mpc_weight_params()
{
    MpcWeightParams w{};

    // ancillary 层相对更关注贴近 z*，同时不过度惩罚控制修正
    const float q[NX] = {2600.0f, 2600.0f, 1400.0f, 0.5f, 0.5f, 0.5f};
    const float r[NU] = {0.08f, 0.08f, 0.08f};
    const float p[NX] = {2600.0f, 2600.0f, 1400.0f, 0.5f, 0.5f, 0.5f};

    for (int i = 0; i < NX; ++i) {
        w.Q[i] = q[i];
        w.P[i] = p[i];
    }

    for (int i = 0; i < NU; ++i) {
        w.R[i] = r[i];
    }

    return w;
}

inline MpcBoundsParams default_central_mpc_bounds_params()
{
    MpcBoundsParams b{};

    // central-path 采用 tightened nominal input set V，为 ancillary 留出调节余量
    b.umin[0] = -28.0f;
    b.umin[1] = -28.0f;
    b.umin[2] =  -7.0f;

    b.umax[0] =  28.0f;
    b.umax[1] =  28.0f;
    b.umax[2] =   7.0f;

    b.pmin[0] = -100.0f;
    b.pmin[1] = -100.0f;
    b.pmin[2] =  -50.0f;

    b.pmax[0] =  100.0f;
    b.pmax[1] =  100.0f;
    b.pmax[2] =   50.0f;

    return b;
}

inline MpcBoundsParams default_ancillary_mpc_bounds_params()
{
    MpcBoundsParams b{};

    // ancillary 层可比 nominal 层更宽，但仍保留到物理限幅之前的一层余量
    b.umin[0] = -40.0f;
    b.umin[1] = -40.0f;
    b.umin[2] =  -9.0f;

    b.umax[0] =  40.0f;
    b.umax[1] =  40.0f;
    b.umax[2] =   9.0f;

    b.pmin[0] = -100.0f;
    b.pmin[1] = -100.0f;
    b.pmin[2] =  -50.0f;

    b.pmax[0] =  100.0f;
    b.pmax[1] =  100.0f;
    b.pmax[2] =   50.0f;

    return b;
}

inline MpcStateBoundsParams default_mpc_state_bounds_params()
{
    MpcStateBoundsParams b{};
    b.enable_soft_bounds = false;

    for (int i = 0; i < NX; ++i) {
        b.xmin[i] = -1.0e6f;
        b.xmax[i] =  1.0e6f;
        b.penalty[i] = 0.0f;
    }

    return b;
}

inline MpcStateBoundsParams default_central_state_bounds_params()
{
    MpcStateBoundsParams b = default_mpc_state_bounds_params();
    b.enable_soft_bounds = true;

    // X/Y 当前不做位置盒约束；先显式约束与安全更直接相关的姿态与速度
    b.penalty[0] = 0.0f;
    b.penalty[1] = 0.0f;

    b.xmin[2] = -2.70f;
    b.xmax[2] =  2.70f;
    b.penalty[2] = 1200.0f;

    b.xmin[3] = -2.40f;
    b.xmax[3] =  2.40f;
    b.penalty[3] = 220.0f;

    b.xmin[4] = -2.40f;
    b.xmax[4] =  2.40f;
    b.penalty[4] = 220.0f;

    b.xmin[5] = -1.60f;
    b.xmax[5] =  1.60f;
    b.penalty[5] = 180.0f;

    return b;
}

inline MpcStateBoundsParams default_ancillary_state_bounds_params()
{
    MpcStateBoundsParams b = default_mpc_state_bounds_params();
    b.enable_soft_bounds = true;

    b.penalty[0] = 0.0f;
    b.penalty[1] = 0.0f;

    b.xmin[2] = -2.95f;
    b.xmax[2] =  2.95f;
    b.penalty[2] = 1800.0f;

    b.xmin[3] = -2.80f;
    b.xmax[3] =  2.80f;
    b.penalty[3] = 320.0f;

    b.xmin[4] = -2.80f;
    b.xmax[4] =  2.80f;
    b.penalty[4] = 320.0f;

    b.xmin[5] = -1.85f;
    b.xmax[5] =  1.85f;
    b.penalty[5] = 260.0f;

    return b;
}

inline MpcConfig default_central_mpc_config()
{
    MpcConfig cfg{};

    cfg.dt = 0.1f;
    cfg.Thor = 1.0f;
    cfg.Nhor = 11;
    cfg.MaxGradIter = 2;
    cfg.MaxMultIter = 1;
    cfg.shift_control = true;
    cfg.use_terminal_cost = true;
    cfg.use_integral_cost = true;

    return cfg;
}

inline MpcConfig default_ancillary_mpc_config()
{
    MpcConfig cfg{};

    // 当前先保持同采样周期/时域长度，后续可独立加快 ancillary 频率。
    cfg.dt = 0.1f;
    cfg.Thor = 1.0f;
    cfg.Nhor = 11;
    cfg.MaxGradIter = 2;
    cfg.MaxMultIter = 1;
    cfg.shift_control = true;
    cfg.use_terminal_cost = true;
    cfg.use_integral_cost = true;

    return cfg;
}

inline MpcBoundsParams paper_reproduction_central_mpc_bounds_params()
{
    MpcBoundsParams b{};

    // 论文给出 nominal tightened bounds = physical - margin
    // 这里对应 [50, 50, 20] - [15, 15, 10] = [35, 35, 10]
    b.umin[0] = -35.0f;
    b.umin[1] = -35.0f;
    b.umin[2] = -10.0f;

    b.umax[0] =  35.0f;
    b.umax[1] =  35.0f;
    b.umax[2] =  10.0f;

    b.pmin[0] = -100.0f;
    b.pmin[1] = -100.0f;
    b.pmin[2] =  -50.0f;

    b.pmax[0] =  100.0f;
    b.pmax[1] =  100.0f;
    b.pmax[2] =   50.0f;

    return b;
}

inline MpcBoundsParams paper_reproduction_ancillary_mpc_bounds_params()
{
    MpcBoundsParams b = paper_reproduction_central_mpc_bounds_params();

    // 当前实现里 ancillary 层仍保留一层接近物理限幅的可调空间；
    // 论文参数尚未把这一层单独展开，因此这里采用保守近物理配置。
    b.umin[0] = -45.0f;
    b.umin[1] = -45.0f;
    b.umin[2] = -18.0f;

    b.umax[0] =  45.0f;
    b.umax[1] =  45.0f;
    b.umax[2] =  18.0f;

    return b;
}

inline RuntimeProfileBundle make_engineering_default_profile_bundle()
{
    RuntimeProfileBundle bundle{};
    bundle.profile_name = "engineering_default";

    bundle.plant = default_plant_params();
    bundle.hopf = default_hopf_config();
    bundle.wave_env = default_wave_env_params();
    bundle.wave_model = default_wave_force_model_params();

    bundle.tube_feedback = default_tube_feedback_params();
    bundle.physical_limits = default_physical_actuator_limits();
    bundle.mission_ref = default_mission_reference_params();
    bundle.experiment = default_experiment_config();

    bundle.central_weights = default_central_mpc_weight_params();
    bundle.central_bounds = default_central_mpc_bounds_params();
    bundle.central_state_bounds = default_central_state_bounds_params();
    bundle.central_config = default_central_mpc_config();

    bundle.ancillary_weights = default_ancillary_mpc_weight_params();
    bundle.ancillary_bounds = default_ancillary_mpc_bounds_params();
    bundle.ancillary_state_bounds = default_ancillary_state_bounds_params();
    bundle.ancillary_config = default_ancillary_mpc_config();

    return bundle;
}

inline RuntimeProfileBundle make_paper_reproduction_profile_bundle()
{
    RuntimeProfileBundle bundle{};
    bundle.profile_name = "paper_reproduction";

    bundle.plant = paper_reproduction_plant_params();
    bundle.hopf = default_hopf_config();
    bundle.wave_env = default_wave_env_params();
    bundle.wave_model = default_wave_force_model_params();

    bundle.tube_feedback = default_tube_feedback_params();
    bundle.physical_limits = paper_reproduction_physical_actuator_limits();
    bundle.mission_ref = paper_reproduction_mission_reference_params();
    bundle.experiment = paper_reproduction_experiment_config();

    bundle.central_weights = default_central_mpc_weight_params();
    bundle.central_bounds = paper_reproduction_central_mpc_bounds_params();
    bundle.central_state_bounds = default_central_state_bounds_params();
    bundle.central_config = default_central_mpc_config();

    bundle.ancillary_weights = default_ancillary_mpc_weight_params();
    bundle.ancillary_bounds = paper_reproduction_ancillary_mpc_bounds_params();
    bundle.ancillary_state_bounds = default_ancillary_state_bounds_params();
    bundle.ancillary_config = default_ancillary_mpc_config();

    return bundle;
}

inline const char *parameter_profile_name(ParameterProfile profile)
{
    switch (profile) {
    case ParameterProfile::PaperReproduction:
        return "paper_reproduction";
    case ParameterProfile::EngineeringDefault:
    default:
        return "engineering_default";
    }
}

inline RuntimeProfileBundle make_runtime_profile_bundle(ParameterProfile profile)
{
    switch (profile) {
    case ParameterProfile::PaperReproduction:
        return make_paper_reproduction_profile_bundle();
    case ParameterProfile::EngineeringDefault:
    default:
        return make_engineering_default_profile_bundle();
    }
}

inline MpcWeightParams default_mpc_weight_params()
{
    MpcWeightParams w{};

    const float q[NX] = {2000.0f, 2000.0f, 1000.0f, 0.2f, 0.2f, 0.2f};
    const float r[NU] = {0.1f, 0.1f, 0.1f};
    const float p[NX] = {2000.0f, 2000.0f, 1000.0f, 0.2f, 0.2f, 0.2f};

    for (int i = 0; i < NX; ++i) {
        w.Q[i] = q[i];
        w.P[i] = p[i];
    }

    for (int i = 0; i < NU; ++i) {
        w.R[i] = r[i];
    }

    return w;
}

inline MpcBoundsParams default_mpc_bounds_params()
{
    MpcBoundsParams b{};

    b.umin[0] = -35.0f;
    b.umin[1] = -35.0f;
    b.umin[2] = -10.0f;

    b.umax[0] =  35.0f;
    b.umax[1] =  35.0f;
    b.umax[2] =  10.0f;

    b.pmin[0] = -100.0f;
    b.pmin[1] = -100.0f;
    b.pmin[2] =  -50.0f;

    b.pmax[0] =  100.0f;
    b.pmax[1] =  100.0f;
    b.pmax[2] =   50.0f;

    return b;
}

inline MpcConfig default_mpc_config()
{
    MpcConfig cfg{};

    cfg.dt = 0.1f;
    cfg.Thor = 1.0f;

    // 与 dt=0.1, Thor=1.0 对应
    cfg.Nhor = 11;

    cfg.MaxGradIter = 2;
    cfg.MaxMultIter = 1;

    cfg.shift_control = true;
    cfg.use_terminal_cost = true;
    cfg.use_integral_cost = true;

    return cfg;
}

} // namespace auv
