#pragma once

#include "auv_common.hpp"

class ActuatorSaturation
{
public:
    using InputVec = auv::InputVec;

    struct DebugInfo {
        InputVec tau_cmd{};
        InputVec tau_applied{};
        InputVec tau_clip{};
    };

    ActuatorSaturation();
    explicit ActuatorSaturation(const auv::PhysicalActuatorLimits &limits);

    void set_limits(const auv::PhysicalActuatorLimits &limits);
    const auv::PhysicalActuatorLimits &limits() const { return _limits; }

    InputVec apply(const InputVec &tau_cmd);
    const DebugInfo &debug() const { return _debug; }

private:
    static float clamp(float val, float lo, float hi);

private:
    auv::PhysicalActuatorLimits _limits{};
    DebugInfo _debug{};
};
