#include "auv_model.hpp"

#include <math.h>
#include <px4_platform_common/defines.h>

/**
 * @brief 默认构造函数
 *
 * 使用统一公共头文件中的默认模型参数。
 */
AuvModel::AuvModel() :
    AuvModel(auv::default_plant_params())
{
}

/**
 * @brief 带参数构造函数
 */
AuvModel::AuvModel(const auv::AuvPlantParams &plant) :
    _plant(plant)
{
}

/**
 * @brief 设置模型参数
 */
void AuvModel::set_plant_params(const auv::AuvPlantParams &plant)
{
    _plant = plant;
}

/**
 * @brief 角度归一化到 [-pi, pi]
 *
 * 这样做的目的是避免航向角在积分过程中持续增长，
 * 导致数值上出现不必要的大角度跳变。
 */
// 第 35 行附近
float AuvModel::wrap_pi(float angle)
{
    float res = fmodf(angle + M_PI_F, 2.0f * M_PI_F);
    if (res < 0.0f) {
        res += 2.0f * M_PI_F;
    }
    return res - M_PI_F;
}

/**
 * @brief 对步长做保护
 *
 * 如果 dt 无效，则回退到 0.02 s（50 Hz）
 * 如果 dt 过大，则截断到 0.2 s
 */
float AuvModel::sanitize_dt(float dt)
{
    if (!PX4_ISFINITE(dt) || dt <= 0.0f) {
        return 0.02f;
    }

    if (dt > 0.2f) {
        return 0.2f;
    }

    return dt;
}

/**
 * @brief 计算 AUV 连续时间状态导数
 *
 * 状态：
 *   x = [X, Y, psi, u, v, r]
 *
 * 输入：
 *   tau = [tau_u, tau_v, tau_r]
 *
 * 扰动：
 *   d_hat = [d_u, d_v, d_r]
 *
 * 模型形式：
 *   eta_dot = R(psi) * nu
 *   nu_dot  = M^{-1} [tau + d_hat - C(nu)nu - D(nu)nu]
 *
 * 这里采用与你现有 AUV 平面 3DOF 模型一致的结构：
 *   surge/sway/yaw 三自由度
 */
AuvModel::StateVec AuvModel::state_derivative(const StateVec &x,
                                              const InputVec &tau,
                                              const InputVec &d_hat) const
{
    StateVec xdot{};
    xdot.zero();

    // -----------------------------
    // 1) 从状态向量中取出各物理量
    // -----------------------------
    const float psi = x(2);
    const float u   = x(3);
    const float v   = x(4);
    const float r   = x(5);

    // -----------------------------
    // 2) 从控制输入中取出广义力
    // -----------------------------
    const float tau_u = tau(0);
    const float tau_v = tau(1);
    const float tau_r = tau(2);

    // -----------------------------
    // 3) 从扰动项中取出各通道补偿/外扰
    // -----------------------------
    const float d_u = d_hat(0);
    const float d_v = d_hat(1);
    const float d_r = d_hat(2);

    // -----------------------------
    // 4) 运动学方程
    //
    // 惯性系位置变化由体坐标速度通过航向角旋转得到：
    //
    //   X_dot = cos(psi) * u - sin(psi) * v
    //   Y_dot = sin(psi) * u + cos(psi) * v
    //   psi_dot = r
    // -----------------------------
    xdot(0) = cosf(psi) * u - sinf(psi) * v;
    xdot(1) = sinf(psi) * u + cosf(psi) * v;
    xdot(2) = r;

    // -----------------------------
    // 5) 动力学方程
    //
    // 这里按照你现在采用的平面 3DOF AUV 结构：
    //
    // surge:
    //   m_u * u_dot = tau_u + d_u - D_u*u + m_v*v*r
    //
    // sway:
    //   m_v * v_dot = tau_v + d_v - D_v*v - m_u*u*r
    //
    // yaw:
    //   Iz * r_dot = tau_r + d_r - D_r*r + (m_u - m_v)*u*v
    //
    // 注意：
    // 这里的符号形式与前面 Hopf 逆动力学重构保持一致。
    // -----------------------------
    xdot(3) = (tau_u + d_u - _plant.Du * u + _plant.m_v * v * r) / _plant.m_u;
    xdot(4) = (tau_v + d_v - _plant.Dv * v - _plant.m_u * u * r) / _plant.m_v;
    xdot(5) = (tau_r + d_r - _plant.Dr * r + (_plant.m_u - _plant.m_v) * u * v) / _plant.Iz;

    return xdot;
}

/**
 * @brief 使用欧拉法推进一步
 *
 * x_{k+1} = x_k + dt * f(x_k, u_k, d_k)
 *
 * 第一版工程实现用欧拉法已经足够清晰和稳定，
 * 后面如果你需要更高精度，可以再升级成 RK4。
 */
AuvModel::StateVec AuvModel::step(const StateVec &x,
                                  const InputVec &tau,
                                  const InputVec &d_hat,
                                  float dt) const
{
    const float h = sanitize_dt(dt);

    // 计算当前导数
    const StateVec xdot = state_derivative(x, tau, d_hat);

    // 欧拉积分
    StateVec x_next = x + h * xdot;

    // 对航向角做归一化，避免角度无限增长
    x_next(2) = wrap_pi(x_next(2));

    return x_next;
}
