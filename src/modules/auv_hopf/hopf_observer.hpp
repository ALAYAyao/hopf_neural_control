#pragma once

#include "auv_common.hpp"

class HopfObserver
{
public:
    using StateVec = auv::StateVec;
    using InputVec = auv::InputVec;

    struct ChannelDebug {
        float d_raw{0.0f};
        float d_pred_before{0.0f};
        float f_error{0.0f};
    };

    struct DebugInfo {
        ChannelDebug surge{};
        ChannelDebug sway{};
        ChannelDebug yaw{};
    };

    struct OscillatorBank {
        float xi[auv::N_OSC]{};
        float yi[auv::N_OSC]{};
        float wi[auv::N_OSC]{};
        float ai[auv::N_OSC]{};
    };

    HopfObserver();
    HopfObserver(const auv::AuvPlantParams &plant, const auv::HopfConfig &hopf_cfg);

    void set_plant_params(const auv::AuvPlantParams &plant);
    void set_hopf_config(const auv::HopfConfig &hopf_cfg);

    InputVec update(const StateVec &x, const InputVec &tau_prev, float Ts);
    bool predict_horizon(float *dhat_hor, int nhor, float Ts) const;

    void reset();

    const InputVec &disturbance_estimate() const { return _d_hat; }
    const DebugInfo &debug() const { return _debug; }
    bool initialized() const { return _initialized; }

private:
    static float safe_sqrt_arg(float v);
    float sanitize_dt(float Ts) const;
    void init_bank(OscillatorBank &bank, const auv::HopfChannelInit &init);
    float reconstruct_disturbance(const OscillatorBank &bank) const;
    void step_bank_euler(OscillatorBank &bank, float F, float Ts, const auv::HopfParams &p);
    void step_bank_preview_autonomous(OscillatorBank &bank, float Ts, const auv::HopfParams &p) const;

private:
    auv::AuvPlantParams _plant{};
    auv::HopfConfig _hopf_cfg{};

    OscillatorBank _osc_u{};
    OscillatorBank _osc_v{};
    OscillatorBank _osc_r{};

    StateVec _x_prev{};
    InputVec _d_hat{};
    bool _initialized{false};
    DebugInfo _debug{};
};
