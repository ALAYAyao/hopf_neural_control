#pragma once

#include "auv_common.hpp"

/**
 * @file eso_observer.hpp
 *
 * @brief 3-DOF AUV velocity-channel linear ESO / LESO disturbance observer.
 *
 * State convention follows auv_common.hpp:
 *   x   = [X, Y, psi, u, v, r]
 *   tau = [tau_u, tau_v, tau_r]
 *   d   = [d_u, d_v, d_r]
 *
 * The observer is converted from the MATLAB LESO block:
 *   e_eso      = z1_eso - nu
 *   nu_dot_est = M^{-1}(F_known + tau_prev + d_hat)
 *   dz1        = nu_dot_est - L1 e_eso
 *   dd_hat     = -L2 e_eso
 *
 * with
 *   L1 = 2 omega_o I,
 *   L2 = omega_o^2 M.
 */
class EsoObserver
{
public:
    using StateVec = auv::StateVec;
    using InputVec = auv::InputVec;

    struct DebugInfo {
        InputVec z1{};          ///< estimated body velocity state [u, v, r]
        InputVec e{};           ///< z1 - measured nu
        InputVec f_known{};     ///< known dynamics term in force/moment coordinates
        InputVec nu_dot_est{};  ///< estimated acceleration of body velocity
        InputVec d_hat_dot{};   ///< disturbance estimate derivative
        bool initialized{false};
    };

    EsoObserver();
    explicit EsoObserver(const auv::AuvPlantParams &plant, float omega_o = 5.0f);

    void set_plant_params(const auv::AuvPlantParams &plant);
    void set_bandwidth(float omega_o);
    float bandwidth() const { return _omega_o; }

    /**
     * @brief Reset observer state.
     *
     * If init_with_measurement is true, the first update initializes z1 to the
     * measured body velocity to avoid a large artificial transient.
     */
    void reset(bool init_with_measurement = true);

    /**
     * @brief One-step ESO update.
     *
     * @param x        current measured AUV state
     * @param tau_prev previously applied actuator force/moment
     * @param Ts       sampling time
     * @return current disturbance estimate d_hat
     */
    InputVec update(const StateVec &x, const InputVec &tau_prev, float Ts);

    /**
     * @brief Fill disturbance horizon by zero-order hold of current d_hat.
     *
     * ESO does not provide a future internal-model rollout, so this method uses
     * the same estimate for all prediction stages.
     */
    bool predict_horizon(float *dhat_hor, int nhor) const;

    const InputVec &disturbance() const { return _d_hat; }
    const DebugInfo &debug() const { return _debug; }

private:
    static float sanitize_dt(float Ts);
    static float safe_mass(float m);
    static float safe_bandwidth(float omega_o);
    static float clampf(float v, float lo, float hi);

private:
    auv::AuvPlantParams _plant{};
    float _omega_o{5.0f};

    InputVec _z1{};
    InputVec _d_hat{};

    bool _initialized{false};
    bool _init_with_measurement{true};

    DebugInfo _debug{};
};
