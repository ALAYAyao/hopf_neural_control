#pragma once

#include "auv_common.hpp"
#include <stdint.h>

class HopfObserver
{
public:
    using StateVec = auv::StateVec;
    using InputVec = auv::InputVec;

    struct ChannelDebug {
        float d_raw;
        float d_pred_before;
        float f_error;
    };

    struct DebugInfo {
        ChannelDebug surge;
        ChannelDebug sway;
        ChannelDebug yaw;
    };

    HopfObserver();
    HopfObserver(const auv::AuvPlantParams &plant, const auv::HopfConfig &hopf_cfg);

    void reset();

    InputVec update(const StateVec &x, const InputVec &tau_prev, float Ts);

    const DebugInfo &debug() const { return _debug; }
    InputVec last_estimate() const { return _d_hat; }

    void set_plant_params(const auv::AuvPlantParams &plant);
    void set_hopf_config(const auv::HopfConfig &hopf_cfg);

private:
    struct OscillatorBank {
        float xi[auv::N_OSC];
        float yi[auv::N_OSC];
        float wi[auv::N_OSC];
        float ai[auv::N_OSC];
    };

    auv::AuvPlantParams _plant{};
    auv::HopfConfig _hopf_cfg{};

    OscillatorBank _osc_u{};
    OscillatorBank _osc_v{};
    OscillatorBank _osc_r{};

    StateVec _x_prev{};
    bool _initialized{false};

    InputVec _d_hat{};
    DebugInfo _debug{};

    void init_bank(OscillatorBank &bank, const auv::HopfChannelInit &init);
    void step_bank_euler(OscillatorBank &bank, float F, float Ts, const auv::HopfParams &p);
    float reconstruct_disturbance(const OscillatorBank &bank) const;

    float sanitize_dt(float Ts) const;
    static float safe_sqrt_arg(float v);
};