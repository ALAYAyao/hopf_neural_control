#pragma once

#include "auv_common.hpp"

class ThrusterAllocation
{
public:
    using StateVec = auv::StateVec;
    using InputVec = auv::InputVec;

    struct DebugInfo {
        InputVec e_vel;
        InputVec pos_error_body;
        InputVec tau_fb;
        InputVec tau_ancillary;
        InputVec tau_raw;
        InputVec tau_applied;
    };

    ThrusterAllocation();
    ThrusterAllocation(const auv::AuvPlantParams &plant,
                       const auv::TubeFeedbackParams &params);

    void set_plant_params(const auv::AuvPlantParams &plant);
    void set_params(const auv::TubeFeedbackParams &params);

    const auv::AuvPlantParams &plant_params() const { return _plant; }
    const auv::TubeFeedbackParams &params() const { return _params; }

    InputVec compute(const StateVec &x,
                     const StateVec &z,
                     const InputVec &d_hat,
                     const InputVec &tau_nom);

    const DebugInfo &debug() const { return _debug; }

private:
    auv::AuvPlantParams _plant{};
    auv::TubeFeedbackParams _params{};
    DebugInfo _debug{};

    static float wrap_pi(float a);
    static float clamp(float val, float lo, float hi);
};