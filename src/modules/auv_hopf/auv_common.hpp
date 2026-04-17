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
static constexpr int NHOR_MAX = 192;

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
 * 当前版本只保留一套工程默认参数。
 * 若后续论文复现需要单独恢复一套口径，再在此基础上扩展即可。
 */
enum class ParameterProfile : uint8_t {
    EngineeringDefault = 0,
};

/**
 * @brief 主模块一次性装载的参数包
 *
 * 当前统一由一套工程默认参数集中装载，避免主模块从多个默认工厂零散拼装，
 * 导致参数口径漂移。
 */

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

struct RuntimeProfileBundle {
    const char *profile_name;

    AuvPlantParams plant;
    HopfConfig hopf;
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

// --- 物理常数 ---
    p.rho = 1025.0f;           // 海水密度 (淡水测试请改为 1000.0f)
    p.g   = 9.81f;             // 重力加速度

    // --- 几何特征 (假设典型的小型框架式 AUV) ---
    p.L     = 0.5f;            // 艇长 50cm
    p.D_auv = 0.3f;            // 等效直径 30cm

    // --- 水动力学特征 ---
    p.cd  = 1.0f;              // 框架式 AUV 阻力较大
    p.cm1 = 0.15f;             // 纵向附加质量系数较小
    p.cm2 = 0.9f;              // 横向附加质量系数较大 (接近圆柱横截面)

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
    p.circle_radius = 30.0f;
    p.circle_omega  = 0.02f;
    return p;
}

inline ExperimentConfig default_experiment_config()
{
    ExperimentConfig cfg{};
    cfg.sim_time_limit = 300.0f;
    return cfg;
}

inline MpcWeightParams default_central_mpc_weight_params()
{
    MpcWeightParams w{};

    const float q[NX] = {800.0f, 800.0f, 500.0f, 0.2f, 0.2f, 0.2f};
    // 增加推力惩罚，阻止满舵飙车
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
    const float q[NX] = {1500.0f, 1500.0f, 1000.0f, 0.5f, 0.5f, 0.5f};
    // 推力惩罚可以低一点，允许在浪大的时候“瞬间发力”
    const float r[NU] = {0.5f, 0.5f, 1.0f};
    const float p[NX] = {1500.0f, 1500.0f, 1000.0f, 0.5f, 0.5f, 0.5f};

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
    b.umin[0] = -45.0f; b.umax[0] =  45.0f;
    b.umin[1] = -35.0f; b.umax[1] =  35.0f;
    b.umin[2] =  -9.0f; b.umax[2] =   9.0f;

    return b;
}

inline MpcBoundsParams default_ancillary_mpc_bounds_params()
{
    MpcBoundsParams b{};

  // --- 1. 控制器输出约束 (umin / umax) ---
    // ancillary 层直接面对真实状态，它的约束需要比 nominal 层（48, 32, 9）更宽。
    // 我们将其设置在非常接近物理极限（60, 40, 12）的位置，留下最后一点防饱和余量。
    b.umin[0] = -55.0f; b.umax[0] =  55.0f;
    b.umin[1] = -45.0f; b.umax[1] =  45.0f;
    b.umin[2] = -11.0f; b.umax[2] =  11.0f;
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
    // 将惩罚设为 0，意味着解除状态空间绝对位置的软约束。
    // (让 MPC 的 Q 矩阵专注于消除与参考点的“相对误差”，而不是限制绝对坐标)
    b.penalty[0] = 0.0f; // X 不限
    b.penalty[1] = 0.0f; // Y 不限
    b.penalty[2] = 0.0f; // psi 不限 (重点修复：允许自由打转和画圆！)

    // --- 2. 纵向速度 (u) ---
    // 小型 AUV 的安全巡航速度通常在 1.0 m/s 左右。
    // 这里限制名义层最高规划速度为 1.2 m/s，确保电机绝对能达标。
    b.xmin[3] = -1.20f;
    b.xmax[3] =  1.20f;
    b.penalty[3] = 250.0f;

    // --- 3. 横向速度 (v) ---
    // 全驱 AUV 虽然能横向平移，但横向迎流面积大，阻力极大。
    // 限制最高横移速度为 0.8 m/s，防止规划出不切实际的侧滑轨迹。
    b.xmin[4] = -0.80f;
    b.xmax[4] =  0.80f;
    b.penalty[4] = 250.0f;

    // --- 4. 艏摇角速度 (r) ---
    // 限制在 1.0 rad/s (约 57度/秒)。
    // 这个速度足够灵活，同时又能防止 AUV 像陀螺一样疯转导致传感器(IMU)晕眩。
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
    // 和名义层一样，绝对不限制全局坐标和航向角，允许 AUV 满世界跑和自由画整圆。
    b.penalty[0] = 0.0f;
    b.penalty[1] = 0.0f;
    b.penalty[2] = 0.0f; // 彻底解除画圆时的 180度/pi 封印

    // --- 2. 纵向速度 (u) 物理极限防线 ---
    // 名义层限速是 1.2m/s。
    // 这里允许遇到大浪时，瞬间被推到 1.5m/s。但超过 1.5m/s 就必须遭受重罚拉回！
    b.xmin[3] = -1.50f;
    b.xmax[3] =  1.50f;
    b.penalty[3] = 300.0f;

    // --- 3. 横向速度 (v) 物理极限防线 ---
    // 名义层规划横移最大 0.8m/s。
    // 真实情况横向受流面积大，允许瞬间漂移到 1.2m/s。
    b.xmin[4] = -1.20f;
    b.xmax[4] =  1.20f;
    b.penalty[4] = 300.0f;

    // --- 4. 艏摇角速度 (r) 物理极限防线 ---
    // 名义层转弯角速度最大 1.0 rad/s。
    // 允许瞬间被海浪拍转到 1.5 rad/s (约85度/秒)，再大就要出安全事故了。
    b.xmin[5] = -1.50f;
    b.xmax[5] =  1.50f;
    b.penalty[5] = 250.0f;

    return b;
}

inline MpcConfig default_central_mpc_config()
{
    MpcConfig cfg{};

    cfg.dt = 0.1f;
    cfg.Thor = 15.0f;
    cfg.Nhor = 151;

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
