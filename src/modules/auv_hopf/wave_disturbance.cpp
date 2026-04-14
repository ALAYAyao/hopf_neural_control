#include "wave_disturbance.hpp"

#include <math.h>
#include <string.h>

using matrix::Vector3f;

WaveDisturbance::WaveDisturbance() :
    WaveDisturbance(auv::default_wave_env_params(),
                    auv::default_wave_force_model_params())
{
}

WaveDisturbance::WaveDisturbance(const auv::WaveEnvParams &wave_env,
                                 const auv::WaveForceModelParams &model_params) :
    _wave_env(wave_env),
    _model_params(model_params)
{
    reset();
}

void WaveDisturbance::reset()
{
    _x_prev_internal.zero();
    _initialized = false;
    _last_force = {};
}

void WaveDisturbance::set_wave_env(const auv::WaveEnvParams &wave_env)
{
    _wave_env = wave_env;
}

void WaveDisturbance::set_model_params(const auv::WaveForceModelParams &model_params)
{
    _model_params = model_params;
}

float WaveDisturbance::wrap_to_pi(float a)
{
    while (a > M_PI_F) {
        a -= 2.0f * M_PI_F;
    }

    while (a < -M_PI_F) {
        a += 2.0f * M_PI_F;
    }

    return a;
}

float WaveDisturbance::trapz(const float *x, const float *y, int n)
{
    if (n <= 1) {
        return 0.0f;
    }

    float s = 0.0f;

    for (int i = 0; i < n - 1; i++) {
        const float dx = x[i + 1] - x[i];
        s += 0.5f * dx * (y[i + 1] + y[i]);
    }

    return s;
}

WaveDisturbance::InputVec WaveDisturbance::step(const StateVec &x, float t)
{
    // 仿真第一步：上一时刻状态初始化为当前状态
    if (!_initialized) {
        _x_prev_internal = x;
        _initialized = true;
    }

    // 取 AUV 速度：来源版只用 [u, v, 0]
    Vector3f U_auv{};
    U_auv(0) = x(3);
    U_auv(1) = x(4);
    U_auv(2) = 0.0f;

    // 来源版里 U_dot_auv 设为零
    Vector3f U_dot_auv{};
    U_dot_auv.zero();

    // 当前遭遇角与上一时刻遭遇角
    const float psi_curr  = x(2);
    const float beta_curr = wrap_to_pi(-psi_curr);
    const float beta_prev = wrap_to_pi(-_x_prev_internal(2));

    // 调用来源版核心波浪力函数
    _last_force = waveforce_sim(_wave_env.omegaMax,
                                _wave_env.Epsilon,
                                _wave_env.Omega,
                                _wave_env.S,
                                _wave_env.mu,
                                U_auv,
                                U_dot_auv,
                                _wave_env.d,
                                beta_prev,
                                beta_curr,
                                t);

    InputVec d_true{};
    d_true(0) = _model_params.scale_factor * _last_force.X;
    d_true(1) = _model_params.scale_factor * _last_force.Y;
    d_true(2) = _model_params.scale_factor * _last_force.N;

    // 更新内部记忆，供下一步使用
    _x_prev_internal = x;

    return d_true;
}

WaveDisturbance::WaveForceOutput WaveDisturbance::waveforce_sim(
    float omegaMax,
    const float Epsilon[auv::N_WAVE_FREQ][auv::N_WAVE_DIR],
    const float Omega[auv::N_WAVE_FREQ],
    const float S[auv::N_WAVE_FREQ],
    const float mu_in[auv::N_WAVE_DIR],
    const Vector3f &U,
    const Vector3f &U_dot,
    float d,
    float beta_last,
    float beta,
    float t) const
{
    WaveForceOutput out{};
    out.X = 0.0f; out.Y = 0.0f; out.Z = 0.0f; out.M = 0.0f; out.N = 0.0f;

    const float rho = _model_params.rho;
    const float L   = _model_params.L;
    const float D   = _model_params.D_auv;
    const float cd  = _model_params.cd;
    const float cm1 = _model_params.cm1;
    const float cm2 = _model_params.cm2;
    const float g   = _model_params.g;
    const float h   = _model_params.h;

    const int Data_num = static_cast<int>(L / h) + 1;
    constexpr int numFreqIntervals = auv::N_WAVE_FREQ;
    constexpr int numDirections    = auv::N_WAVE_DIR;

    const float deltaOmega = omegaMax / static_cast<float>(numFreqIntervals);
    float deltaDirections = 1.0f;
    if (numDirections != 1) {
        deltaDirections = static_cast<float>(M_PI) / static_cast<float>(numDirections - 1);
    }

    // 1) 处理方向角 mu
    float mu[auv::N_WAVE_DIR];
    for (int i = 0; i < numDirections; i++) {
        mu[i] = mu_in[i] - beta_last + beta;
    }

    // 2) 方向扩散函数 M(mu)
    float Mdir[auv::N_WAVE_DIR];
    if (numDirections != 1) {
        for (int i = 0; i < numDirections; i++) {
            Mdir[i] = (2.0f / static_cast<float>(M_PI)) * cosf(mu[i]) * cosf(mu[i]);
        }
        const float integral_M = trapz(mu, Mdir, numDirections);
        if (fabsf(integral_M) > 1e-8f) {
            for (int i = 0; i < numDirections; i++) {
                Mdir[i] /= integral_M;
            }
        }
    } else {
        Mdir[0] = 1.0f;
    }

    // 3) 计算振幅 Amp
    float Amp[auv::N_WAVE_FREQ][auv::N_WAVE_DIR];
    for (int n = 0; n < numFreqIntervals; n++) {
        for (int m = 0; m < numDirections; m++) {
            const float val = 2.0f * (S[n] * Mdir[m]) * deltaOmega * deltaDirections;
            Amp[n][m] = (val > 0.0f) ? sqrtf(val) : 0.0f;
        }
    }

    // 4) 核心物理面积常数
    const float Ax = static_cast<float>(M_PI) * D * D / 4.0f; // 截面积
    const float Ay = Ax; // 假设横截面积相同

    // ==============================================================
    // 5) 核心合并循环：去除局部大数组，在同一个切片内边算速度、边算受力
    // ==============================================================
    for (int j = 0; j < Data_num; j++) {
        const float x_pos = static_cast<float>(j) * h - L / 2.0f;
        const float y_pos = 0.0f;

        // --- A. 计算当前切片流体速度 (替代原 v_horizontal 等数组) ---
        float v_horiz = 0.0f, v_horiz_dot = 0.0f;
        float w_vert = 0.0f,  w_vert_dot = 0.0f;

        for (int n = 0; n < numFreqIntervals; n++) {
            for (int m = 0; m < numDirections; m++) {
                const float omega_enc = Omega[n] - (Omega[n] * Omega[n] / g) * U(0) * cosf(mu[m]);
                const float phase = t * omega_enc - (Omega[n] * Omega[n] / g) * (x_pos * cosf(mu[m]) + y_pos * sinf(mu[m])) + Epsilon[n][m];
                const float atten = expf(-(Omega[n] * Omega[n] / g) * d);

                v_horiz      +=  Omega[n] * atten * Amp[n][m] * cosf(phase);
                v_horiz_dot  -=  Omega[n] * Omega[n] * atten * Amp[n][m] * sinf(phase);
                w_vert       +=  Omega[n] * atten * Amp[n][m] * sinf(phase);
                w_vert_dot   +=  Omega[n] * Omega[n] * atten * Amp[n][m] * cosf(phase);
            }
        }

        // --- B. 速度分解到船体坐标系 ---
        const float u_seg = v_horiz * cosf(beta);
        const float v_seg = v_horiz * sinf(beta);
        const float u_dot_seg = v_horiz_dot * cosf(beta);
        const float v_dot_seg = v_horiz_dot * sinf(beta);

        // --- C. 相对速度与相对加速度 ---
        const float u_rel = u_seg - U(0);
        const float v_rel = v_seg - U(1);
        const float w_rel = w_vert - U(2);

        const float u_rel_dot = u_dot_seg - U_dot(0);
        const float v_rel_dot = v_dot_seg - U_dot(1);
        const float w_rel_dot = w_vert_dot - U_dot(2);

        // --- D. 计算切片受力 (严格还原 MATLAB 逻辑) ---
        
        // 1. X 轴受力 (附加质量项每个切片都有)
        float dFx = rho * cm1 * Ax * u_rel_dot * h;
        // 严格还原 MATLAB: 阻力只加在首尾两端 (j==0 或 j==Data_num-1)
        if (j == 0 || j == Data_num - 1) {
            dFx += 0.5f * rho * cd * Ax * (u_rel * fabsf(u_rel) + U(0) * fabsf(U(0)));
        }

        // 2. Y 轴受力 (附加质量 + 圆柱侧向阻力)
        const float dFy = rho * cm2 * Ay * v_rel_dot * h 
                        + 0.5f * rho * cd * D * v_rel * fabsf(v_rel) * h;

        // 3. Z 轴受力
        const float dFz = rho * cm2 * Ay * w_rel_dot * h 
                        + 0.5f * rho * cd * D * w_rel * fabsf(w_rel) * h;

        // --- E. 积分累加 ---
        out.X += dFx;
        out.Y += dFy;
        out.Z += dFz;
        out.M += x_pos * dFz; // 纵摇力矩
        out.N += x_pos * dFy; // 艏摇力矩
    }

    return out;
}