#pragma once

/**
 * @file auv_model.hpp
 *
 * @brief AUV 名义动力学模型
 *
 * 该类负责：
 * 1. 保存 AUV 的动力学参数
 * 2. 根据当前状态 x 和控制输入 tau 计算状态导数 x_dot
 * 3. 使用简单欧拉积分执行一步离散推进
 *
 * 状态定义：
 *   x = [X, Y, psi, u, v, r]
 *
 * 控制定义：
 *   tau = [tau_u, tau_v, tau_r]
 *
 * 其中：
 *   X, Y   : 惯性坐标系位置
 *   psi    : 航向角
 *   u, v   : 体坐标系纵向/横向速度
 *   r      : 偏航角速度
 */

#include "auv_common.hpp"

class AuvModel
{
public:
    using StateVec = auv::StateVec;
    using InputVec = auv::InputVec;

    /**
     * @brief 默认构造，使用 auv_common.hpp 中的默认模型参数
     */
    AuvModel();

    /**
     * @brief 使用指定模型参数构造
     * @param plant AUV 动力学参数
     */
    explicit AuvModel(const auv::AuvPlantParams &plant);

    /**
     * @brief 设置模型参数
     * @param plant 新的 AUV 动力学参数
     */
    void set_plant_params(const auv::AuvPlantParams &plant);

    /**
     * @brief 获取当前模型参数
     */
    const auv::AuvPlantParams &plant_params() const { return _plant; }

    /**
     * @brief 计算连续时间状态导数
     *
     * @param x     当前状态 [X, Y, psi, u, v, r]
     * @param tau   当前控制输入 [tau_u, tau_v, tau_r]
     * @param d_hat 外部扰动估计/补偿项
     *
     * @return x_dot = f(x, tau, d_hat)
     *
     * 说明：
     * 这里把 d_hat 作为“外部等效广义力”加入动力学。
     * 如果你想计算纯名义模型，可以传入零向量。
     */
    StateVec state_derivative(const StateVec &x,
                              const InputVec &tau,
                              const InputVec &d_hat) const;

    /**
     * @brief 使用欧拉法做一步离散推进
     *
     * @param x     当前状态
     * @param tau   控制输入
     * @param d_hat 扰动项
     * @param dt    步长（秒）
     *
     * @return 下一时刻状态
     */
    StateVec step(const StateVec &x,
                  const InputVec &tau,
                  const InputVec &d_hat,
                  float dt) const;

private:
    /// AUV 动力学参数
    auv::AuvPlantParams _plant{};

    /**
     * @brief 将角度归一化到 [-pi, pi]
     */
    static float wrap_pi(float angle);

    /**
     * @brief 对积分步长做保护，避免异常 dt 导致数值爆炸
     */
    static float sanitize_dt(float dt);
};