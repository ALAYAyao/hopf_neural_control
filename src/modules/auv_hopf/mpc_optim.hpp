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
        double uref_hor[auv::NU * auv::NHOR_MAX]{};
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

    /**
     * @brief 设置带状态参考、控制参考和扰动预览的完整时域
     *
     * 对 central-path OCP：
     *   xref_hor = mission reference
     *   uref_hor = nullptr / zero
     *   dhat_hor = zero preview
     *
     * 对 ancillary NMPC：
     *   xref_hor = z* horizon
     *   uref_hor = v* horizon
     *   dhat_hor = d_hat[k+i|k]
     */
    void set_reference_horizon_control_preview(const StateVec &z,
                                               const float *xref_hor,
                                               const float *uref_hor,
                                               const float *dhat_hor,
                                               int nhor);

    InputVec solve();
    void reset();

    bool initialized() const { return _initialized; }
    const DebugInfo &debug() const { return _debug; }

    int last_solution_horizon_length() const { return _last_solution_nhor; }

    /**
     * @brief 导出最近一次求解得到的状态/控制时域
     *
     * x_hor 和 u_hor 均按 stage 平铺，长度由 nhor_max 指定。
     * 若 nhor_max 大于实际缓存长度，则尾部用最后一个 stage 保持。
     */
    bool export_last_solution_horizon(float *x_hor,
                                      float *u_hor,
                                      int nhor_max) const;

private:
    void free_grampc();
    void sync_userparam_from_common();
    int compute_nhor(double dt, double Thor) const;
    void cache_last_solution_horizon();

private:
    typeGRAMPC *_grampc{nullptr};

    auv::AuvPlantParams _plant{};
    auv::MpcWeightParams _weights{};
    auv::MpcBoundsParams _bounds{};
    auv::MpcConfig _config{};

    UserParam _userparam{};
    bool _initialized{false};
    DebugInfo _debug{};

    float _last_x_hor[auv::NX * auv::NHOR_MAX]{};
    float _last_u_hor[auv::NU * auv::NHOR_MAX]{};
    int _last_solution_nhor{0};
};
