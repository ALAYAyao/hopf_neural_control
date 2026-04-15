#include "actuator_saturation.hpp"

ActuatorSaturation::ActuatorSaturation() :
    ActuatorSaturation(auv::default_physical_actuator_limits())
{
}

ActuatorSaturation::ActuatorSaturation(const auv::PhysicalActuatorLimits &limits) :
    _limits(limits)
{
    _debug = {};
}

void ActuatorSaturation::set_limits(const auv::PhysicalActuatorLimits &limits)
{
    _limits = limits;
}

float ActuatorSaturation::clamp(float val, float lo, float hi)
{
    if (val > hi) {
        return hi;
    }

    if (val < lo) {
        return lo;
    }

    return val;
}

ActuatorSaturation::InputVec ActuatorSaturation::apply(const InputVec &tau_cmd)
{
    InputVec tau_applied{};

    tau_applied(0) = clamp(tau_cmd(0), _limits.tau_u_min, _limits.tau_u_max);
    tau_applied(1) = clamp(tau_cmd(1), _limits.tau_v_min, _limits.tau_v_max);
    tau_applied(2) = clamp(tau_cmd(2), _limits.tau_r_min, _limits.tau_r_max);

    _debug.tau_cmd = tau_cmd;
    _debug.tau_applied = tau_applied;
    _debug.tau_clip(0) = tau_applied(0) - tau_cmd(0);
    _debug.tau_clip(1) = tau_applied(1) - tau_cmd(1);
    _debug.tau_clip(2) = tau_applied(2) - tau_cmd(2);

    return tau_applied;
}
