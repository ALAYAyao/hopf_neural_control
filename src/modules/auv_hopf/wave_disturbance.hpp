#pragma once

/**
 * @file wave_disturbance.hpp
 *
 * @brief 来源版海浪扰动模块
 *
 * 该类对应原始 MATLAB 结构：
 *
 *   d_true = wave_disturbance(x, t)
 *
 * 内部调用：
 *
 *   waveforce_SIM(...)
 *
 * 特点：
 * 1. 保留 wave_env 参数结构
 * 2. 保留上一时刻状态记忆，用于 beta_last
 * 3. 保留沿艇体长度离散积分
 * 4. 输出真实扰动 d_true = [X_force, Y_force, N_moment]
 *
 * 注意：
 * 这个类输出的是“真实环境外扰”
 * 与 HopfObserver 的 d_hat 不同。
 */

#include "auv_common.hpp"

class WaveDisturbance
{
public:
    using StateVec = auv::StateVec;
    using InputVec = auv::InputVec;

    /**
     * @brief 波浪载荷核心输出
     */
    struct WaveForceOutput {
        float X;
        float Y;
        float Z;
        float M;
        float N;
    };

    WaveDisturbance();
    WaveDisturbance(const auv::WaveEnvParams &wave_env,
                    const auv::WaveForceModelParams &model_params);

    /**
     * @brief 重置内部状态记忆
     */
    void reset();

    /**
     * @brief 设置波浪环境参数
     */
    void set_wave_env(const auv::WaveEnvParams &wave_env);

    /**
     * @brief 设置波浪力模型参数
     */
    void set_model_params(const auv::WaveForceModelParams &model_params);

    /**
     * @brief 获取当前波浪环境参数
     */
    const auv::WaveEnvParams &wave_env() const { return _wave_env; }

    /**
     * @brief 获取当前波浪力模型参数
     */
    const auv::WaveForceModelParams &model_params() const { return _model_params; }

    /**
     * @brief 来源版外层接口
     *
     * @param x 当前状态 [X, Y, psi, u, v, r]
     * @param t 当前时间
     * @return d_true = [X_force, Y_force, N_moment]
     */
    InputVec step(const StateVec &x, float t);

    /**
     * @brief 获取上一次完整波浪力输出（便于调试）
     */
    const WaveForceOutput &last_force_output() const { return _last_force; }

private:
    float _last_time{0.0f};
    float _wave_phase[auv::N_WAVE_FREQ][auv::N_WAVE_DIR]{};
    auv::WaveEnvParams _wave_env{};
    auv::WaveForceModelParams _model_params{};

    StateVec _x_prev_internal{};
    bool _initialized{false};

    WaveForceOutput _last_force{};

    /**
     * @brief 来源版核心：计算波浪力和力矩
     */
    WaveForceOutput waveforce_sim(float omegaMax,
                                  const float Epsilon[auv::N_WAVE_FREQ][auv::N_WAVE_DIR],
                                  const float Omega[auv::N_WAVE_FREQ],
                                  const float S[auv::N_WAVE_FREQ],
                                  const float mu_in[auv::N_WAVE_DIR],
                                  const matrix::Vector3f &U,
                                  const matrix::Vector3f &U_dot,
                                  float d,
                                  float beta_last,
                                  float beta,
                                  float dt);

    /**
     * @brief 角度归一化到 [-pi, pi]
     */
    static float wrap_to_pi(float a);

    /**
     * @brief 梯形积分（等价 MATLAB trapz(x,y)）
     */
    static float trapz(const float *x, const float *y, int n);
};