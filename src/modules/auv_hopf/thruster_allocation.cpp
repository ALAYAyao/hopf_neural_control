#include "thruster_allocation.hpp"

#include <math.h>

ThrusterAllocation::ThrusterAllocation() :
    ThrusterAllocation(auv::default_plant_params(),
                       auv::default_tube_feedback_params())
{
}

ThrusterAllocation::ThrusterAllocation(const auv::AuvPlantParams &plant,
                                       const auv::TubeFeedbackParams &params) :
    _plant(plant),
    _params(params)
{
    _debug = {};
}

void ThrusterAllocation::set_plant_params(const auv::AuvPlantParams &plant)
{
    _plant = plant;
}

void ThrusterAllocation::set_params(const auv::TubeFeedbackParams &params)
{
    _params = params;
}

float ThrusterAllocation::wrap_pi(float a)
{
// 对应 MATLAB 的 mod(a + pi, 2*pi) - pi
    float wrapped = fmodf(a + M_PI_F, 2.0f * M_PI_F);
    if (wrapped < 0.0f) {
        wrapped += 2.0f * M_PI_F; // fmodf 对负数求余可能仍是负数，修正一下
    }
    return wrapped - M_PI_F;
}

float ThrusterAllocation::clamp(float val, float lo, float hi)
{
    if (val > hi) {
        return hi;
    }

    if (val < lo) {
        return lo;
    }

    return val;
}

ThrusterAllocation::InputVec ThrusterAllocation::compute(const StateVec &x,
                                                         const StateVec &z,
                                                         const InputVec &d_hat,
                                                         const InputVec &tau_nom)
{
    // --------------------------------------------------
    // 1) 误差 e = x - z
    // --------------------------------------------------
    StateVec e = x - z;
    e(2) = wrap_pi(e(2));

    // --------------------------------------------------
    // 2) 质量/惯量矩阵 M（来自 AuvPlantParams）
    //    M = diag([m_u, m_v, Iz])
    // --------------------------------------------------
    const float M11 = _plant.m_u;
    const float M22 = _plant.m_v;
    const float M33 = _plant.Iz;

    // --------------------------------------------------
    // 3) Kp, Kv（来自 TubeFeedbackParams）
    // --------------------------------------------------
    const float Kp11 = _params.Kp11;
    const float Kp22 = _params.Kp22;
    const float Kp33 = _params.Kp33;

    const float Kv11 = _params.Kv11;
    const float Kv22 = _params.Kv22;
    const float Kv33 = _params.Kv33;

    // --------------------------------------------------
    // 4) tau_fb = -M * (Kv * e(4:6) + Kp * [0;0;e(3)])
    // --------------------------------------------------
    InputVec tau_fb{};
    tau_fb(0) = -M11 * (Kv11 * e(3) + Kp11 * 0.0f);
    tau_fb(1) = -M22 * (Kv22 * e(4) + Kp22 * 0.0f);
    tau_fb(2) = -M33 * (Kv33 * e(5) + Kp33 * e(2));

    // --------------------------------------------------
    // 5) 位置误差转体坐标
    // --------------------------------------------------
    const float psi = x(2);

    InputVec pos_error_body{};
    pos_error_body(0) =  cosf(psi) * e(0) + sinf(psi) * e(1);
    pos_error_body(1) = -sinf(psi) * e(0) + cosf(psi) * e(1);
    pos_error_body(2) = 0.0f;

    // --------------------------------------------------
    // 6) tau_ancillary = tau_fb - M * Kp * pos_error_body
    // --------------------------------------------------
    InputVec tau_ancillary{};
    tau_ancillary(0) = tau_fb(0) - M11 * Kp11 * pos_error_body(0);
    tau_ancillary(1) = tau_fb(1) - M22 * Kp22 * pos_error_body(1);
    tau_ancillary(2) = tau_fb(2) - M33 * Kp33 * pos_error_body(2);

    // --------------------------------------------------
    // 7) 总控制
    // --------------------------------------------------
    InputVec tau_raw{};
    tau_raw(0) = tau_nom(0) + tau_ancillary(0) - d_hat(0);
    tau_raw(1) = tau_nom(1) + tau_ancillary(1) - d_hat(1);
    tau_raw(2) = tau_nom(2) + tau_ancillary(2) - d_hat(2);

    // --------------------------------------------------
    // 8) 执行层物理限幅
    // --------------------------------------------------
    InputVec tau_applied{};
    tau_applied(0) = clamp(tau_raw(0), -_params.tau_u_max_phys, _params.tau_u_max_phys);
    tau_applied(1) = clamp(tau_raw(1), -_params.tau_v_max_phys, _params.tau_v_max_phys);
    tau_applied(2) = clamp(tau_raw(2), -_params.tau_r_max_phys, _params.tau_r_max_phys);

    // --------------------------------------------------
    // 9) 调试量
    // --------------------------------------------------
    _debug.e_vel(0) = e(3);
    _debug.e_vel(1) = e(4);
    _debug.e_vel(2) = e(5);

    _debug.pos_error_body = pos_error_body;
    _debug.tau_fb = tau_fb;
    _debug.tau_ancillary = tau_ancillary;
    _debug.tau_raw = tau_raw;
    _debug.tau_applied = tau_applied;

    return tau_applied;
}
