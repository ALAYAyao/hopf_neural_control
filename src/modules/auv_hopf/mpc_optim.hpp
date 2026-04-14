#pragma once

#include "auv_common.hpp"

extern "C" {
#include "grampc.h"
}

class MpcOptim
{
public:
    using StateVec = auv::StateVec;
    using InputVec = auv::InputVec;
    using ParamVec = auv::ParamVec;

    struct DebugInfo {
        bool initialized{false};
        float dt{0.0f};
        float Thor{0.0f};
        int Nhor{0};
        float t0{0.0f};
        InputVec last_u_nom{};
    };

    // ======================================================
    // 公开给 mpc_optim.cpp 内部 GRAMPC C 回调桥接使用
    // 这样就不需要再绕 private 可见性问题
    // ======================================================
    struct UserParam {
        double m_u{0.0};
        double m_v{0.0};
        double Iz{0.0};

        double Du{0.0};
        double Dv{0.0};
        double Dr{0.0};

        double Q[auv::NX]{};
        double R[auv::NU]{};
        double P[auv::NX]{};

        double xref_hor[auv::NX * auv::NHOR_MAX]{};
        double dhat_hor[auv::NP * auv::NHOR_MAX]{};
        int Nhor_ref{2};
    };

    MpcOptim();
    MpcOptim(const auv::AuvPlantParams &plant,
             const auv::MpcWeightParams &weights,
             const auv::MpcBoundsParams &bounds,
             const auv::MpcConfig &config);

    int horizon_length() const { return _debug.Nhor; }
    float sample_time() const { return _debug.dt; }
    float horizon_time() const { return _debug.Thor; }

    void initialize(float dt, float Thor);

    void set_plant_params(const auv::AuvPlantParams &plant);
    void set_weight_params(const auv::MpcWeightParams &weights);
    void set_bounds_params(const auv::MpcBoundsParams &bounds);
    void set_mpc_config(const auv::MpcConfig &config);

    const auv::AuvPlantParams &plant_params() const { return _plant; }
    const auv::MpcWeightParams &weight_params() const { return _weights; }
    const auv::MpcBoundsParams &bounds_params() const { return _bounds; }
    const auv::MpcConfig &config() const { return _config; }

    void set_reference_horizon(const StateVec &z,
                               const float *xref_hor,
                               int nhor,
                               const ParamVec &dhat);

    /**
     * @brief 设置带多步扰动预览的参考时域
     *
     * xref_hor 与 dhat_hor 都按 stage 平铺：
     *   xref_hor = [x0(0..NX-1), x1(...), ...]
     *   dhat_hor = [d0(0..NP-1), d1(...), ...]
     *
     * 这是当前从“单步 d_hat 占位”过渡到“预测层真实 preview”所需的关键接口。
     */
    void set_reference_horizon_preview(const StateVec &z,
                                       const float *xref_hor,
                                       const float *dhat_hor,
                                       int nhor);

    InputVec solve();
    void reset();

    bool initialized() const { return _initialized; }
    const DebugInfo &debug() const { return _debug; }

private:
    void free_grampc();
    void sync_userparam_from_common();
    int compute_nhor(double dt, double Thor) const;

private:
    typeGRAMPC *_grampc{nullptr};

    auv::AuvPlantParams _plant{};
    auv::MpcWeightParams _weights{};
    auv::MpcBoundsParams _bounds{};
    auv::MpcConfig _config{};

    UserParam _userparam{};
    bool _initialized{false};
    DebugInfo _debug{};
};