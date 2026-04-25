#pragma once

/**
 * @file auv_common.hpp
 *
 * @brief AUV Hopf / ESO + Tube-MPC 项目的公共类型、参数与默认配置
 *
 * 本文件是整个项目的统一公共入口，负责集中定义：
 * 1. 基本维度与公共向量类型
 * 2. AUV 名义动力学参数
 * 3. Hopf 学习器参数
 * 4. ESO 观测器参数
 * 5. 波浪环境与波浪力模型参数
 * 6. MPC / GRAMPC 参数、约束与默认配置
 *
 * 说明：
 * - 为了保证 Hopf 与 ESO 对比实验公平，plant / wave / actuator / MPC 参数保持完全一致。
 * - ESO 只额外增加自身带宽 omega_o，不改变原工程参数口径。
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
// ESO 观测器参数
// ==================================================

/**
 * @brief ESO 对比实验参数
 *
 * 仅用于配置 LESO 带宽，不改变 plant / wave / MPC / actuator 参数。
 * 这样 Hopf 与 ESO 对比实验只替换 disturbance observer / preview mechanism。
 */
struct EsoConfig {
    float omega_o;
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
// 执行层 / 任务层参数
// ==================================================

/**
 * @brief 执行层物理限幅参数
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
 * @brief 运行实验配置
 */
struct ExperimentConfig {
    float sim_time_limit;
};

/**
 * @brief 参数口径分组
 *
 * 当前版本只保留一套工程默认参数。
 * 若后续论文复现需要单独恢复一套口径，再在此基础上扩展即可。
 */
enum class ParameterProfile : uint8_t {
    EngineeringDefault = 0,
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

struct MissionReferenceParams {
    // 从离线参考轨迹表的哪个点开始
    int ref_start_index;

    // 每个控制拍推进几个参考点，通常为 1
    int ref_stride;
};

/**
 * @brief 主模块一次性装载的参数包
 *
 * 当前统一由一套工程默认参数集中装载，避免主模块从多个默认工厂零散拼装，
 * 导致参数口径漂移。
 */
struct RuntimeProfileBundle {
    const char *profile_name;

    AuvPlantParams plant;
    HopfConfig hopf;
    EsoConfig eso;
    WaveEnvParams wave_env;
    WaveForceModelParams wave_model;

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
// 默认工厂函数
// ==================================================

inline AuvPlantParams default_plant_params()
{
    AuvPlantParams p{};

    // 质量与附加质量 (单位: kg)
    p.m_u = 25.0f;  // 纵向质量 = 20kg(干重) + 5kg(附加质量)
    p.m_v = 35.0f;  // 横向质量 = 20kg(干重) + 15kg(附加质量)
    p.Iz  = 3.5f;   // 转动惯量 (约 0.5米长的方形/圆柱体，含附加惯量)

    // 线性阻尼 (单位: N/(m/s) 和 Nm/(rad/s))
    p.Du  = 12.0f;  // 纵向流体阻尼
    p.Dv  = 25.0f;  // 横向流体阻尼 (横向面积更大，阻尼更大)
    p.Dr  = 8.0f;   // 艏摇旋转阻尼

    return p;
}

inline HopfConfig default_hopf_config()
{
    HopfConfig cfg{};

    cfg.surge.gamma = 1.0f;
    cfg.surge.mu = 1.0f;
    cfg.surge.epsF = 1.7f;
    cfg.surge.epsW = 0.55f;
    cfg.surge.eta = 1.0f;
    cfg.surge.lambda_alpha = 0.0f;

    cfg.sway.gamma = 1.0f;
    cfg.sway.mu = 1.0f;
    cfg.sway.epsF = 1.7f;
    cfg.sway.epsW = 0.55f;
    cfg.sway.eta = 0.5f;
    cfg.sway.lambda_alpha = 0.0f;

    // --- Yaw 通道 (低通抗噪模式) ---
    cfg.yaw.gamma = 1.0f;
    cfg.yaw.mu = 1.0f;

    // 因为有 LPF 护体，直接恢复和 Surge/Sway 一样的满血学习率。
    cfg.yaw.epsF = 1.7f;
    cfg.yaw.epsW = 0.55f;
    cfg.yaw.eta = 1.5f;
    cfg.yaw.lambda_alpha = 0.0f;

    const float w0[N_OSC]     = {0.5f, 1.3333f, 2.1667f, 3.0f};
    const float w_yaw[N_OSC] = {0.5f, 1.3333f, 0.0f, 0.0f};

    const float xu0[N_OSC] = { 0.0291f,  0.0083f,  0.0138f, -0.0106f};
    const float yu0[N_OSC] = {-0.0047f, -0.0027f,  0.0110f, -0.0028f};

    const float xv0[N_OSC] = { 0.0070f, -0.0205f, -0.0035f, -0.0082f};
    const float yv0[N_OSC] = {-0.0158f,  0.0051f,  0.0028f,  0.0003f};

    const float xr0[N_OSC] = {-3.00f,  2.50f,  0.0000f,  0.0000f};
    const float yr0[N_OSC] = { 1.00f, -1.00f,  0.0000f,  0.0000f};

    for (int i = 0; i < N_OSC; ++i) {
        cfg.init_surge.xi0[i] = xu0[i];
        cfg.init_surge.yi0[i] = yu0[i];
        cfg.init_surge.wi0[i] = w0[i];

        cfg.init_sway.xi0[i] = xv0[i];
        cfg.init_sway.yi0[i] = yv0[i];
        cfg.init_sway.wi0[i] = w0[i];

        cfg.init_yaw.xi0[i] = xr0[i];
        cfg.init_yaw.yi0[i] = yr0[i];
        cfg.init_yaw.wi0[i] = w_yaw[i];
    }

    return cfg;
}

inline EsoConfig default_eso_config()
{
    EsoConfig cfg{};

    // 与 MATLAB LESO 形式一致的默认带宽；只影响 ESO，对 Hopf 工程参数无影响。
    cfg.omega_o = 5.0f;

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

    // --- 物理常数 ---
    p.rho = 1025.0f;           // 海水密度 (淡水测试请改为 1000.0f)
    p.g   = 9.81f;             // 重力加速度

    // --- 几何特征 (假设典型的小型框架式 AUV) ---
    p.L     = 0.5f;            // 艇长 50cm
    p.D_auv = 0.3f;            // 等效直径 30cm

    // --- 水动力学特征 ---
    p.cd  = 1.0f;              // 框架式 AUV 阻力较大
    p.cm1 = 0.15f;             // 纵向附加质量系数较小
    p.cm2 = 0.9f;              // 横向附加惯性相关系数较大 (接近圆柱横截面)

    // --- 算法参数 ---
    p.h = p.L / 10.0f;         // 离散步长 0.05m (分成 10 个切片)
    p.scale_factor = 1.0f;     // 初始波浪扰动缩放保持原样

    return p;
}

inline PhysicalActuatorLimits default_physical_actuator_limits()
{
    PhysicalActuatorLimits p{};

    p.tau_u_min = -60.0f; p.tau_u_max =  60.0f;
    p.tau_v_min = -50.0f; p.tau_v_max =  50.0f;
    p.tau_r_min = -12.0f; p.tau_r_max =  12.0f;

    return p;
}

inline MissionReferenceParams default_mission_reference_params()
{
    MissionReferenceParams p{};
    p.ref_start_index = 0;
    p.ref_stride = 1;
    return p;
}

inline ExperimentConfig default_experiment_config()
{
    ExperimentConfig cfg{};
    cfg.sim_time_limit = 150.0f;
    return cfg;
}

inline MpcWeightParams default_central_mpc_weight_params()
{
    MpcWeightParams w{};

    const float q[NX] = {1200.0f, 1200.0f, 1200.0f, 0.2f, 0.2f, 0.2f};
    const float r[NU] = {2.0f, 2.0f, 5.0f};
    const float p[NX] = {800.0f, 800.0f, 500.0f, 0.2f, 0.2f, 0.2f};

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
    const float q[NX] = {1200.0f, 1200.0f, 1200.0f, 0.2f, 0.2f, 0.2f};
    const float r[NU] = {2.0f, 2.0f, 5.0f};
    const float p[NX] = {800.0f, 800.0f, 500.0f, 0.2f, 0.2f, 0.2f};

    for (int i = 0; i < NX; ++i) {
        w.Q[i] = q[i];
        w.P[i] = p[i];
    }

    for (int i = 0; i < NU; ++i) {
        w.R[i] = r[i];
    }

    return w;
}

inline void set_default_disturbance_parameter_bounds(MpcBoundsParams &b)
{
    // GRAMPC 的 p = [d_u, d_v, d_r] 作为外生扰动估计序列传入。
    // 当前代码主要通过 dhat_hor 使用扰动，不主动优化 p。
    // 这里显式给出宽边界，避免未初始化边界影响求解器。
    b.pmin[0] = -500.0f; b.pmax[0] = 500.0f;
    b.pmin[1] = -500.0f; b.pmax[1] = 500.0f;
    b.pmin[2] = -500.0f; b.pmax[2] = 500.0f;
}

inline MpcBoundsParams default_central_mpc_bounds_params()
{
    MpcBoundsParams b{};

    // central-path 采用 tightened nominal input set V，为 ancillary 留出调节余量
    b.umin[0] = -45.0f; b.umax[0] =  45.0f;
    b.umin[1] = -35.0f; b.umax[1] =  35.0f;
    b.umin[2] =  -9.0f; b.umax[2] =   9.0f;

    set_default_disturbance_parameter_bounds(b);

    return b;
}

inline MpcBoundsParams default_ancillary_mpc_bounds_params()
{
    MpcBoundsParams b{};

    // ancillary 层直接面对真实状态，它的约束需要比 nominal 层更宽。
    // 设置在接近物理极限的位置，留下最后一点防饱和余量。
    b.umin[0] = -55.0f; b.umax[0] =  55.0f;
    b.umin[1] = -45.0f; b.umax[1] =  45.0f;
    b.umin[2] = -11.0f; b.umax[2] =  11.0f;

    set_default_disturbance_parameter_bounds(b);

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

    // --- 1. 绝对坐标与姿态 (X, Y, psi) ---
    // AUV 需要全海域航行（X/Y无界），并且要能画整圆（psi无界）。
    b.penalty[0] = 0.0f;
    b.penalty[1] = 0.0f;
    b.penalty[2] = 0.0f;

    // --- 2. 纵向速度 (u) ---
    b.xmin[3] = -1.20f;
    b.xmax[3] =  1.20f;
    b.penalty[3] = 250.0f;

    // --- 3. 横向速度 (v) ---
    b.xmin[4] = -0.80f;
    b.xmax[4] =  0.80f;
    b.penalty[4] = 250.0f;

    // --- 4. 艏摇角速度 (r) ---
    b.xmin[5] = -1.00f;
    b.xmax[5] =  1.00f;
    b.penalty[5] = 200.0f;

    return b;
}

inline MpcStateBoundsParams default_ancillary_state_bounds_params()
{
    MpcStateBoundsParams b = default_mpc_state_bounds_params();
    b.enable_soft_bounds = true;

    // --- 1. 位置与姿态 (X, Y, psi) ---
    b.penalty[0] = 0.0f;
    b.penalty[1] = 0.0f;
    b.penalty[2] = 0.0f;

    // --- 2. 纵向速度 (u) 物理极限防线 ---
    b.xmin[3] = -1.50f;
    b.xmax[3] =  1.50f;
    b.penalty[3] = 300.0f;

    // --- 3. 横向速度 (v) 物理极限防线 ---
    b.xmin[4] = -1.20f;
    b.xmax[4] =  1.20f;
    b.penalty[4] = 300.0f;

    // --- 4. 艏摇角速度 (r) 物理极限防线 ---
    b.xmin[5] = -1.50f;
    b.xmax[5] =  1.50f;
    b.penalty[5] = 250.0f;

    return b;
}

inline MpcConfig default_central_mpc_config()
{
    MpcConfig cfg{};

    cfg.dt = 0.1f;
    cfg.Thor = 1.5f;
    cfg.Nhor = 16;

    cfg.MaxGradIter = 10;
    cfg.MaxMultIter = 3;
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
    cfg.MaxGradIter = 10;
    cfg.MaxMultIter = 3;
    cfg.shift_control = true;
    cfg.use_terminal_cost = true;
    cfg.use_integral_cost = true;

    return cfg;
}

inline RuntimeProfileBundle make_runtime_profile_bundle()
{
    RuntimeProfileBundle bundle{};
    bundle.profile_name = "engineering_default";

    bundle.plant = default_plant_params();
    bundle.hopf = default_hopf_config();
    bundle.eso = default_eso_config();
    bundle.wave_env = default_wave_env_params();
    bundle.wave_model = default_wave_force_model_params();

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

// 兼容旧调用接口：当前只有一套工程默认参数，因此 profile 参数被忽略。
inline const char *parameter_profile_name(ParameterProfile)
{
    return "engineering_default";
}

inline RuntimeProfileBundle make_runtime_profile_bundle(ParameterProfile)
{
    return make_runtime_profile_bundle();
}

inline MpcWeightParams default_mpc_weight_params()
{
    return default_ancillary_mpc_weight_params();
}

inline MpcBoundsParams default_mpc_bounds_params()
{
    return default_ancillary_mpc_bounds_params();
}

inline MpcConfig default_mpc_config()
{
    return default_ancillary_mpc_config();
}

} // namespace auv
