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

    /**
     * @brief 基于当前已学得的 Hopf 内模状态，生成短时扰动预览序列
     *
     * 语义说明：
     * 1. 这里生成的是 learned internal model 的 short-horizon signal extrapolation；
     * 2. 它不是基于未来海浪真参数或未来状态-海浪耦合关系的真扰动预测；
     * 3. 因而该预览应被理解为 ancillary MPC 的外生 preview 输入，
     *    未被该内模捕获的误差部分仍由 residual / tube 机制处理。
     */
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
    // 仅用于对已学得振子内模做短时自主前滚，不代表真实海浪机理的未来精确传播。
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
