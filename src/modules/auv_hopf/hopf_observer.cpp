#include "hopf_observer.hpp"
#include <math.h>
#include <px4_platform_common/defines.h>

HopfObserver::HopfObserver() :
    HopfObserver(auv::default_plant_params(), auv::default_hopf_config())
{
}

HopfObserver::HopfObserver(const auv::AuvPlantParams &plant, const auv::HopfConfig &hopf_cfg) :
    _plant(plant),
    _hopf_cfg(hopf_cfg)
{
    reset();
}

float HopfObserver::safe_sqrt_arg(float v)
{
    return (v < 1e-8f) ? 1e-8f : v;
}

float HopfObserver::sanitize_dt(float Ts) const
{
    if (!PX4_ISFINITE(Ts) || Ts <= 0.0f) {
        return 0.02f;
    }

    if (Ts > 0.2f) {
        return 0.2f;
    }

    return Ts;
}

void HopfObserver::set_plant_params(const auv::AuvPlantParams &plant)
{
    _plant = plant;
}

void HopfObserver::set_hopf_config(const auv::HopfConfig &hopf_cfg)
{
    _hopf_cfg = hopf_cfg;
    reset();
}

void HopfObserver::init_bank(OscillatorBank &bank, const auv::HopfChannelInit &init)
{
    for (int i = 0; i < auv::N_OSC; i++) {
        bank.xi[i] = init.xi0[i];
        bank.yi[i] = init.yi0[i];
        bank.wi[i] = init.wi0[i];
        bank.ai[i] = 0.0f;
    }
}

void HopfObserver::reset()
{
    init_bank(_osc_u, _hopf_cfg.init_surge);
    init_bank(_osc_v, _hopf_cfg.init_sway);
    init_bank(_osc_r, _hopf_cfg.init_yaw);

    _x_prev.zero();
    _d_hat.zero();
    _initialized = false;
    _debug = {};
}

float HopfObserver::reconstruct_disturbance(const OscillatorBank &bank) const
{
    float sum = 0.0f;

    for (int i = 0; i < auv::N_OSC; i++) {
        sum += bank.ai[i] * bank.xi[i];
    }

    return sum;
}

// 该步骤只是在冻结当前已学习内模的前提下，对振子组做短时自主前滚。
// 它表达的是 learned internal model 的外推，而不是未来真海浪-真遭遇频率的精确机理预测。
void HopfObserver::step_bank_preview_autonomous(OscillatorBank &bank, float Ts, const auv::HopfParams &p) const
{
    for (int i = 0; i < auv::N_OSC; i++) {
        const float x = bank.xi[i];
        const float y = bank.yi[i];
        const float w = bank.wi[i];
        const float a = bank.ai[i];

        const float r2 = x * x + y * y;

        const float dx = p.gamma * (p.mu - r2) * x - w * y;
        const float dy = p.gamma * (p.mu - r2) * y + w * x;
        const float dw = 0.0f;
        const float da = -p.lambda_alpha * a;

        bank.xi[i] += Ts * dx;
        bank.yi[i] += Ts * dy;
        bank.wi[i] += Ts * dw;
        bank.ai[i] += Ts * da;

        if (!PX4_ISFINITE(bank.xi[i])) { bank.xi[i] = 0.0f; }
        if (!PX4_ISFINITE(bank.yi[i])) { bank.yi[i] = 0.0f; }
        if (!PX4_ISFINITE(bank.wi[i])) { bank.wi[i] = 0.5f; }
        if (!PX4_ISFINITE(bank.ai[i])) { bank.ai[i] = 0.0f; }

        if (bank.wi[i] < 0.1f) {
            bank.wi[i] = 0.1f;
        }

        if (bank.wi[i] > 5.0f) {
            bank.wi[i] = 5.0f;
        }

        if (bank.ai[i] > 10.0f) {
            bank.ai[i] = 10.0f;
        } else if (bank.ai[i] < -10.0f) {
            bank.ai[i] = -10.0f;
        }
    }
}

// 生成的 dhat_hor 用于 ancillary MPC 的外生扰动预览输入。
// 这里复制当前振子 bank 并做短时前滚，得到的是“当前已学得扰动信号内模”的外推序列。
// 因为控制器并不知道未来海况真值，所以这里并不尝试构造状态耦合的真扰动预报。
bool HopfObserver::predict_horizon(float *dhat_hor, int nhor, float Ts) const
{
    if (dhat_hor == nullptr) {
        return false;
    }

    const int nh = (nhor < 2) ? 2 : ((nhor > auv::NHOR_MAX) ? auv::NHOR_MAX : nhor);
    const float dt = sanitize_dt(Ts);

    for (int k = 0; k < auv::NHOR_MAX; ++k) {
        for (int i = 0; i < auv::NP; ++i) {
            dhat_hor[k * auv::NP + i] = 0.0f;
        }
    }

    if (!_initialized) {
        return false;
    }

    OscillatorBank osc_u = _osc_u;
    OscillatorBank osc_v = _osc_v;
    OscillatorBank osc_r = _osc_r;

    for (int k = 0; k < nh; ++k) {
        dhat_hor[k * auv::NP + 0] = reconstruct_disturbance(osc_u);
        dhat_hor[k * auv::NP + 1] = reconstruct_disturbance(osc_v);
        dhat_hor[k * auv::NP + 2] = reconstruct_disturbance(osc_r);

        step_bank_preview_autonomous(osc_u, dt, _hopf_cfg.surge);
        step_bank_preview_autonomous(osc_v, dt, _hopf_cfg.sway);
        step_bank_preview_autonomous(osc_r, dt, _hopf_cfg.yaw);
    }

    for (int k = nh; k < auv::NHOR_MAX; ++k) {
        dhat_hor[k * auv::NP + 0] = dhat_hor[(nh - 1) * auv::NP + 0];
        dhat_hor[k * auv::NP + 1] = dhat_hor[(nh - 1) * auv::NP + 1];
        dhat_hor[k * auv::NP + 2] = dhat_hor[(nh - 1) * auv::NP + 2];
    }

    return true;
}

void HopfObserver::step_bank_euler(OscillatorBank &bank, float F, float Ts, const auv::HopfParams &p)
{
    for (int i = 0; i < auv::N_OSC; i++) {
        const float x = bank.xi[i];
        const float y = bank.yi[i];
        const float w = bank.wi[i];
        const float a = bank.ai[i];

        const float r2 = x * x + y * y;
        const float r2_safe = safe_sqrt_arg(r2);

        const float dx = p.gamma * (p.mu - r2) * x - w * y + p.epsF * F;
        const float dy = p.gamma * (p.mu - r2) * y + w * x;
        const float dw = -p.epsW * F * (y / sqrtf(r2_safe));
        const float da = p.eta * x * F - p.lambda_alpha * a;

        bank.xi[i] += Ts * dx;
        bank.yi[i] += Ts * dy;
        bank.wi[i] += Ts * dw;
        bank.ai[i] += Ts * da;

        if (!PX4_ISFINITE(bank.xi[i])) { bank.xi[i] = 0.0f; }
        if (!PX4_ISFINITE(bank.yi[i])) { bank.yi[i] = 0.0f; }
        if (!PX4_ISFINITE(bank.wi[i])) { bank.wi[i] = 0.5f; }
        if (!PX4_ISFINITE(bank.ai[i])) { bank.ai[i] = 0.0f; }

        if (bank.wi[i] < 0.1f) {
            bank.wi[i] = 0.1f;
        }

        if (bank.wi[i] > 20.0f) {
            bank.wi[i] = 20.0f;
        }

        if (bank.ai[i] > 200.0f) {
            bank.ai[i] = 200.0f;
        } else if (bank.ai[i] < -200.0f) {
            bank.ai[i] = -200.0f;
        }
    }
}

HopfObserver::InputVec HopfObserver::update(const StateVec &x, const InputVec &tau_prev, float Ts)
{
    const float dt = sanitize_dt(Ts);

    if (!_initialized) {
        _x_prev = x;
        _initialized = true;
        _d_hat.zero();
        _debug = {};
        return _d_hat;
    }

    const float u_prev = _x_prev(3);
    const float v_prev = _x_prev(4);
    const float r_prev = _x_prev(5);

    const float du = (x(3) - _x_prev(3)) / dt;
    const float dv = (x(4) - _x_prev(4)) / dt;
    const float dr = (x(5) - _x_prev(5)) / dt;

    const float F_surge = _plant.m_u * du + _plant.Du * u_prev - _plant.m_v * v_prev * r_prev;
    const float F_sway  = _plant.m_v * dv + _plant.Dv * v_prev + _plant.m_u * u_prev * r_prev;
    const float F_yaw   = _plant.Iz  * dr + _plant.Dr * r_prev - (_plant.m_u - _plant.m_v) * u_prev * v_prev;

    InputVec d_raw;
    d_raw(0) = F_surge - tau_prev(0);
    d_raw(1) = F_sway  - tau_prev(1);
    d_raw(2) = F_yaw   - tau_prev(2); // 因为 yaw 通道的扰动被大幅衰减了，所以在喂给 Hopf 前先放大 100 倍，让它更容易学习和捕获 yaw 方向的扰动特征。
    // static float d_raw_yaw_filtered = 0.0f;
    // if (!_initialized) {
    //     d_raw_yaw_filtered = d_raw(2);
    // }
    // // 滤波系数 alpha = 0.1，只信任 10% 的新数据，保留 90% 的历史平滑数据
    // d_raw_yaw_filtered = 0.9f * d_raw_yaw_filtered + 0.1f * d_raw(2);

    // // 把滤波后的干净信号喂给 Hopf
    // d_raw(2) = d_raw_yaw_filtered;

    const float d_pred_u = reconstruct_disturbance(_osc_u);
    const float d_pred_v = reconstruct_disturbance(_osc_v);
    const float d_pred_r = reconstruct_disturbance(_osc_r);

    const float F_err_u = d_raw(0) - d_pred_u;
    const float F_err_v = d_raw(1) - d_pred_v;
    const float F_err_r = d_raw(2) - d_pred_r;

    _debug.surge.d_raw = d_raw(0);
    _debug.surge.d_pred_before = d_pred_u;
    _debug.surge.f_error = F_err_u;

    _debug.sway.d_raw = d_raw(1);
    _debug.sway.d_pred_before = d_pred_v;
    _debug.sway.f_error = F_err_v;

    _debug.yaw.d_raw = d_raw(2);
    _debug.yaw.d_pred_before = d_pred_r;
    _debug.yaw.f_error = F_err_r;

    step_bank_euler(_osc_u, F_err_u, dt, _hopf_cfg.surge);
    step_bank_euler(_osc_v, F_err_v, dt, _hopf_cfg.sway);
    step_bank_euler(_osc_r, F_err_r, dt, _hopf_cfg.yaw);

    _d_hat(0) = reconstruct_disturbance(_osc_u);
    _d_hat(1) = reconstruct_disturbance(_osc_v);
    _d_hat(2) = reconstruct_disturbance(_osc_r);

    for (int i = 0; i < auv::NU; i++) {
        if (!PX4_ISFINITE(_d_hat(i))) {
            _d_hat(i) = 0.0f;
        }
    }

    _x_prev = x;
    return _d_hat;
}
