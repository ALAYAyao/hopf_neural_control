#include "eso_observer.hpp"

#include <math.h>
#include <px4_platform_common/defines.h>

EsoObserver::EsoObserver() :
    EsoObserver(auv::default_plant_params(), 5.0f)
{
}

EsoObserver::EsoObserver(const auv::AuvPlantParams &plant, float omega_o) :
    _plant(plant),
    _omega_o(safe_bandwidth(omega_o))
{
    reset(true);
}

void EsoObserver::set_plant_params(const auv::AuvPlantParams &plant)
{
    _plant = plant;
}

void EsoObserver::set_bandwidth(float omega_o)
{
    _omega_o = safe_bandwidth(omega_o);
}

void EsoObserver::reset(bool init_with_measurement)
{
    _z1.zero();
    _d_hat.zero();
    _initialized = false;
    _init_with_measurement = init_with_measurement;
    _debug = {};
}

float EsoObserver::sanitize_dt(float Ts)
{
    if (!PX4_ISFINITE(Ts) || Ts <= 0.0f) {
        return 0.02f;
    }

    if (Ts > 0.2f) {
        return 0.2f;
    }

    return Ts;
}

float EsoObserver::safe_mass(float m)
{
    return (fabsf(m) < 1e-6f) ? 1.0f : m;
}

float EsoObserver::safe_bandwidth(float omega_o)
{
    if (!PX4_ISFINITE(omega_o) || omega_o <= 0.0f) {
        return 5.0f;
    }

    // 防止离散欧拉更新下过高带宽导致数值发散。
    if (omega_o > 50.0f) {
        return 50.0f;
    }

    return omega_o;
}

float EsoObserver::clampf(float v, float lo, float hi)
{
    if (v < lo) { return lo; }
    if (v > hi) { return hi; }
    return v;
}

EsoObserver::InputVec EsoObserver::update(const StateVec &x,
                                          const InputVec &tau_prev,
                                          float Ts)
{
    const float dt = sanitize_dt(Ts);

    const float u = x(3);
    const float v = x(4);
    const float r = x(5);

    InputVec nu{};
    nu(0) = u;
    nu(1) = v;
    nu(2) = r;

    if (!_initialized) {
        if (_init_with_measurement) {
            _z1 = nu;
        }
        _initialized = true;
    }

    const float m_u = safe_mass(_plant.m_u);
    const float m_v = safe_mass(_plant.m_v);
    const float Iz  = safe_mass(_plant.Iz);

    // MATLAB 对应：
    // F_known = [-Du*u + m_v*v*r;
    //            -Dv*v - m_u*u*r;
    //            -Dr*r + (m_u - m_v)*u*v]
    InputVec f_known{};
    f_known(0) = -_plant.Du * u + _plant.m_v * v * r;
    f_known(1) = -_plant.Dv * v - _plant.m_u * u * r;
    f_known(2) = -_plant.Dr * r + (_plant.m_u - _plant.m_v) * u * v;

    InputVec e{};
    e(0) = _z1(0) - nu(0);
    e(1) = _z1(1) - nu(1);
    e(2) = _z1(2) - nu(2);

    InputVec nu_dot_est{};
    nu_dot_est(0) = (f_known(0) + tau_prev(0) + _d_hat(0)) / m_u;
    nu_dot_est(1) = (f_known(1) + tau_prev(1) + _d_hat(1)) / m_v;
    nu_dot_est(2) = (f_known(2) + tau_prev(2) + _d_hat(2)) / Iz;

    const float L1 = 2.0f * _omega_o;
    const float wo2 = _omega_o * _omega_o;

    InputVec dz1{};
    dz1(0) = nu_dot_est(0) - L1 * e(0);
    dz1(1) = nu_dot_est(1) - L1 * e(1);
    dz1(2) = nu_dot_est(2) - L1 * e(2);

    InputVec d_hat_dot{};
    d_hat_dot(0) = -wo2 * m_u * e(0);
    d_hat_dot(1) = -wo2 * m_v * e(1);
    d_hat_dot(2) = -wo2 * Iz  * e(2);

    _z1(0) += dt * dz1(0);
    _z1(1) += dt * dz1(1);
    _z1(2) += dt * dz1(2);

    _d_hat(0) += dt * d_hat_dot(0);
    _d_hat(1) += dt * d_hat_dot(1);
    _d_hat(2) += dt * d_hat_dot(2);

    // 工程保护：避免异常数值污染优化器。
    for (int i = 0; i < auv::NU; ++i) {
        if (!PX4_ISFINITE(_z1(i))) {
            _z1(i) = nu(i);
        }
        if (!PX4_ISFINITE(_d_hat(i))) {
            _d_hat(i) = 0.0f;
        }
        _d_hat(i) = clampf(_d_hat(i), -500.0f, 500.0f);
    }

    _debug.z1 = _z1;
    _debug.e = e;
    _debug.f_known = f_known;
    _debug.nu_dot_est = nu_dot_est;
    _debug.d_hat_dot = d_hat_dot;
    _debug.initialized = _initialized;

    return _d_hat;
}

bool EsoObserver::predict_horizon(float *dhat_hor, int nhor) const
{
    if (dhat_hor == nullptr) {
        return false;
    }

    int nh = nhor;
    if (nh < 2) {
        nh = 2;
    }
    if (nh > auv::NHOR_MAX) {
        nh = auv::NHOR_MAX;
    }

    for (int k = 0; k < auv::NHOR_MAX; ++k) {
        const bool active = (k < nh);
        for (int i = 0; i < auv::NP; ++i) {
            dhat_hor[k * auv::NP + i] = active ? _d_hat(i) : _d_hat(i);
        }
    }

    return _initialized;
}
