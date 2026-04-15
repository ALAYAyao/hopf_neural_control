#include "mpc_optim.hpp"

#include <math.h>
#include <string.h>

namespace
{

// ======================================================
// 当前活动实例桥接到 GRAMPC C 风格回调
//
// 注意：
// 1) 这里仍然通过全局 g_up 把 C 风格回调绑定到当前实例；
// 2) 但 solve / set_reference... 会在调用前显式激活所属实例，
//    因而可以串行复用多个 MpcOptim（例如 central + ancillary）。
// ======================================================
static MpcOptim::UserParam *g_up = nullptr;

static int clamp_int(int x, int xmin, int xmax)
{
    if (x < xmin) { return xmin; }
    if (x > xmax) { return xmax; }
    return x;
}

static double wrap_pi(double a)
{
    while (a > M_PI)  { a -= 2.0 * M_PI; }
    while (a < -M_PI) { a += 2.0 * M_PI; }
    return a;
}

static int get_stage_index(double t,
                           const typeGRAMPCparam *param,
                           const MpcOptim::UserParam *up)
{
    const double tau = (t - param->t0) / param->dt;
    const int k = static_cast<int>(tau + 0.5);
    return clamp_int(k, 0, up->Nhor_ref - 1);
}

static const double *get_xref_stage(double t,
                                    const typeGRAMPCparam *param,
                                    const MpcOptim::UserParam *up)
{
    const int k = get_stage_index(t, param, up);
    return &up->xref_hor[auv::NX * k];
}

static const double *get_xref_terminal(const MpcOptim::UserParam *up)
{
    const int k = clamp_int(up->Nhor_ref - 1, 0, auv::NHOR_MAX - 1);
    return &up->xref_hor[auv::NX * k];
}


static const double *get_uref_stage(double t,
                                    const typeGRAMPCparam *param,
                                    const MpcOptim::UserParam *up)
{
    const int k = get_stage_index(t, param, up);
    return &up->uref_hor[auv::NU * k];
}

static const double *get_dhat_stage(double t,
                                    const typeGRAMPCparam *param,
                                    const MpcOptim::UserParam *up)
{
    const int k = get_stage_index(t, param, up);
    return &up->dhat_hor[auv::NP * k];
}

} // anonymous namespace

// ======================================================
// GRAMPC required problem functions
//
// 当前实现对应：discrete-time system
// 不使用 semi-implicit ODE / DAE / RODAS 路线，
// 因此不再提供 Mfct / Mtrans 等附加函数。
// ======================================================

extern "C"
{

void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh,
             typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
    (void)userparam;
    *Nx  = auv::NX;
    *Nu  = auv::NU;
    *Np  = auv::NP;
    *Ng  = 0;
    *Nh  = 0;
    *NgT = 0;
    *NhT = 0;
}

void ffct(typeRNum *out,
          ctypeRNum t,
          ctypeRNum *x,
          ctypeRNum *u,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)t;
    (void)userparam;

    const auto *up = g_up;
    const double Ts = param->dt;

    const double psi = x[2];
    const double uu  = x[3];
    const double vv  = x[4];
    const double rr  = x[5];

    const double tau_u = u[0];
    const double tau_v = u[1];
    const double tau_r = u[2];

    const double *dhat = get_dhat_stage(t, param, up);

    const double d_u = dhat[0];
    const double d_v = dhat[1];
    const double d_r = dhat[2];

    const double cpsi = cos(psi);
    const double spsi = sin(psi);

    const double Xdot = cpsi * uu - spsi * vv;
    const double Ydot = spsi * uu + cpsi * vv;
    const double psidot = rr;

    const double udot =
        (tau_u + d_u - up->Du * uu + up->m_v * vv * rr) / up->m_u;

    const double vdot =
        (tau_v + d_v - up->Dv * vv - up->m_u * uu * rr) / up->m_v;

    const double rdot =
        (tau_r + d_r - up->Dr * rr + (up->m_u - up->m_v) * uu * vv) / up->Iz;

    // discrete-time step map
    out[0] = x[0] + Ts * Xdot;
    out[1] = x[1] + Ts * Ydot;
    out[2] = wrap_pi(x[2] + Ts * psidot);
    out[3] = x[3] + Ts * udot;
    out[4] = x[4] + Ts * vdot;
    out[5] = x[5] + Ts * rdot;
}

void dfdx_vec(typeRNum *out,
              ctypeRNum t,
              ctypeRNum *x,
              ctypeRNum *u,
              ctypeRNum *p,
              ctypeRNum *vec,
              const typeGRAMPCparam *param,
              typeUSERPARAM *userparam)
{
    (void)t;
    (void)u;
    (void)p;
    (void)userparam;

    const auto *up = g_up;
    const double Ts = param->dt;

    const double psi = x[2];
    const double uu  = x[3];
    const double vv  = x[4];
    const double rr  = x[5];

    const double cpsi = cos(psi);
    const double spsi = sin(psi);

    // Jacobian of discrete map x_{k+1} = f(x_k,u_k,p_k)
    double A[auv::NX][auv::NX] = {};

    for (int i = 0; i < auv::NX; ++i) {
        A[i][i] = 1.0;
    }

    A[0][2] += Ts * (-spsi * uu - cpsi * vv);
    A[0][3] += Ts * cpsi;
    A[0][4] += Ts * (-spsi);

    A[1][2] += Ts * ( cpsi * uu - spsi * vv);
    A[1][3] += Ts * spsi;
    A[1][4] += Ts * cpsi;

    A[2][5] += Ts;

    A[3][3] += Ts * (-up->Du / up->m_u);
    A[3][4] += Ts * ( up->m_v * rr / up->m_u);
    A[3][5] += Ts * ( up->m_v * vv / up->m_u);

    A[4][3] += Ts * (-up->m_u * rr / up->m_v);
    A[4][4] += Ts * (-up->Dv / up->m_v);
    A[4][5] += Ts * (-up->m_u * uu / up->m_v);

    A[5][3] += Ts * ((up->m_u - up->m_v) * vv / up->Iz);
    A[5][4] += Ts * ((up->m_u - up->m_v) * uu / up->Iz);
    A[5][5] += Ts * (-up->Dr / up->Iz);

    for (int j = 0; j < auv::NX; ++j) {
        out[j] = 0.0;
        for (int i = 0; i < auv::NX; ++i) {
            out[j] += A[i][j] * vec[i];
        }
    }
}

void dfdu_vec(typeRNum *out,
              ctypeRNum t,
              ctypeRNum *x,
              ctypeRNum *u,
              ctypeRNum *p,
              ctypeRNum *vec,
              const typeGRAMPCparam *param,
              typeUSERPARAM *userparam)
{
    (void)t;
    (void)x;
    (void)u;
    (void)p;
    (void)userparam;

    const auto *up = g_up;
    const double Ts = param->dt;

    // B matrix of discrete map
    // x_{k+1} = ... + B u_k
    double B[auv::NX][auv::NU] = {};

    B[3][0] = Ts / up->m_u;
    B[4][1] = Ts / up->m_v;
    B[5][2] = Ts / up->Iz;

    for (int j = 0; j < auv::NU; ++j) {
        out[j] = 0.0;
        for (int i = 0; i < auv::NX; ++i) {
            out[j] += B[i][j] * vec[i];
        }
    }
}

void dfdp_vec(typeRNum *out,
              ctypeRNum t,
              ctypeRNum *x,
              ctypeRNum *u,
              ctypeRNum *p,
              ctypeRNum *vec,
              const typeGRAMPCparam *param,
              typeUSERPARAM *userparam)
{
    (void)t;
    (void)u;
    (void)p;
    (void)userparam;

    const auto *up = g_up;
    const double Ts = param->dt;

    // parameter p = [d_u, d_v, d_r]
    double G[auv::NX][auv::NP] = {};

    G[3][0] = Ts / up->m_u;
    G[4][1] = Ts / up->m_v;
    G[5][2] = Ts / up->Iz;

    for (int j = 0; j < auv::NP; ++j) {
        out[j] = 0.0;
        for (int i = 0; i < auv::NX; ++i) {
            out[j] += G[i][j] * vec[i];
        }
    }
}

void lfct(typeRNum *out,
          ctypeRNum t,
          ctypeRNum *x,
          ctypeRNum *u,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)p;
    (void)userparam;

    const auto *up = g_up;
    const double *xdes = get_xref_stage(t, param, up);

    double cost = 0.0;

    for (int i = 0; i < auv::NX; ++i) {
        double e = x[i] - xdes[i];
        if (i == 2) {
            e = wrap_pi(e);
        }
        cost += 0.5 * up->Q[i] * e * e;
    }

    const double *udes = get_uref_stage(t, param, up);

    for (int i = 0; i < auv::NU; ++i) {
        const double du = u[i] - udes[i];
        cost += 0.5 * up->R[i] * du * du;
    }

    *out = cost;
}

void dldx(typeRNum *out,
          ctypeRNum t,
          ctypeRNum *x,
          ctypeRNum *u,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)u;
    (void)p;
    (void)userparam;

    const auto *up = g_up;
    const double *xdes = get_xref_stage(t, param, up);

    for (int i = 0; i < auv::NX; ++i) {
        double e = x[i] - xdes[i];
        if (i == 2) {
            e = wrap_pi(e);
        }
        out[i] = up->Q[i] * e;
    }
}

void dldu(typeRNum *out,
          ctypeRNum t,
          ctypeRNum *x,
          ctypeRNum *u,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)t;
    (void)x;
    (void)p;
    (void)param;
    (void)userparam;

    const auto *up = g_up;
    const double *udes = get_uref_stage(t, param, up);

    for (int i = 0; i < auv::NU; ++i) {
        out[i] = up->R[i] * (u[i] - udes[i]);
    }
}

void dldp(typeRNum *out,
          ctypeRNum t,
          ctypeRNum *x,
          ctypeRNum *u,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)t;
    (void)x;
    (void)u;
    (void)p;
    (void)param;
    (void)userparam;

    for (int i = 0; i < auv::NP; ++i) {
        out[i] = 0.0;
    }
}

void Vfct(typeRNum *out,
          ctypeRNum T,
          ctypeRNum *x,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)T;
    (void)p;
    (void)param;
    (void)userparam;

    const auto *up = g_up;
    const double *xdes = get_xref_terminal(up);

    double cost = 0.0;
    for (int i = 0; i < auv::NX; ++i) {
        double e = x[i] - xdes[i];
        if (i == 2) {
            e = wrap_pi(e);
        }
        cost += 0.5 * up->P[i] * e * e;
    }

    *out = cost;
}

void dVdx(typeRNum *out,
          ctypeRNum T,
          ctypeRNum *x,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)T;
    (void)p;
    (void)param;
    (void)userparam;

    const auto *up = g_up;
    const double *xdes = get_xref_terminal(up);

    for (int i = 0; i < auv::NX; ++i) {
        double e = x[i] - xdes[i];
        if (i == 2) {
            e = wrap_pi(e);
        }
        out[i] = up->P[i] * e;
    }
}

void dVdp(typeRNum *out,
          ctypeRNum T,
          ctypeRNum *x,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)T;
    (void)x;
    (void)p;
    (void)param;
    (void)userparam;

    for (int i = 0; i < auv::NP; ++i) {
        out[i] = 0.0;
    }
}

void dVdT(typeRNum *out,
          ctypeRNum T,
          ctypeRNum *x,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)T;
    (void)x;
    (void)p;
    (void)param;
    (void)userparam;
    out[0] = 0.0;
}

// ======================================================
// no-constraint stubs
// even when Ng=Nh=NgT=NhT=0, GRAMPC still expects
// these global symbols to exist at link time.
// ======================================================

void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p,
          const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)t; (void)x; (void)u; (void)p; (void)param; (void)userparam;
}

void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p,
          const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)t; (void)x; (void)u; (void)p; (void)param; (void)userparam;
}

void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p,
           const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)T; (void)x; (void)p; (void)param; (void)userparam;
}

void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p,
           const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)T; (void)x; (void)p; (void)param; (void)userparam;
}

void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p,
              ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)t; (void)x; (void)u; (void)p; (void)vec; (void)param; (void)userparam;
}

void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p,
              ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)t; (void)x; (void)u; (void)p; (void)vec; (void)param; (void)userparam;
}

void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p,
              ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)t; (void)x; (void)u; (void)p; (void)vec; (void)param; (void)userparam;
}

void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p,
              ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)t; (void)x; (void)u; (void)p; (void)vec; (void)param; (void)userparam;
}

void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p,
              ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)t; (void)x; (void)u; (void)p; (void)vec; (void)param; (void)userparam;
}

void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p,
              ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)t; (void)x; (void)u; (void)p; (void)vec; (void)param; (void)userparam;
}

void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p,
               ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)T; (void)x; (void)p; (void)vec; (void)param; (void)userparam;
}

void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p,
               ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)T; (void)x; (void)p; (void)vec; (void)param; (void)userparam;
}

void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p,
               ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)x; (void)p; (void)param; (void)userparam;
    out[0] = 0.0;
    (void)T; (void)vec;
}

void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p,
               ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)T; (void)x; (void)p; (void)vec; (void)param; (void)userparam;
}

void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p,
               ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)out; (void)T; (void)x; (void)p; (void)vec; (void)param; (void)userparam;
}

void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p,
               ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)x; (void)p; (void)param; (void)userparam;
    out[0] = 0.0;
    (void)T; (void)vec;
}

// ======================================================
// RODAS compatibility stubs
//
// 说明：
// 1) 当前 AUV MPC 实际采用 discrete-time GRAMPC 配置；
// 2) 这些函数不是当前控制语义所必需；
// 3) 但只要 grampc 库把 rodas.c 链进来，链接器就要求这些符号存在。
// 因此这里提供最保守、最安全的占位实现：
// - 质量矩阵取单位阵
// - 时间导数取 0
// - Jacobian 相关项取 0
//
// 若未来真的切到 semi-implicit ODE / DAE + RODAS，再把这里替换成
// 真实实现。
// ======================================================

void Mfct(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    (void)param;
    (void)userparam;

    // full identity mass matrix, column-wise
    for (int j = 0; j < auv::NX; ++j) {
        for (int i = 0; i < auv::NX; ++i) {
            out[j * auv::NX + i] = (i == j) ? 1.0 : 0.0;
        }
    }
}

void Mtrans(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
{
    Mfct(out, param, userparam);
}

void dfdx(typeRNum *out,
          ctypeRNum t,
          ctypeRNum *x,
          ctypeRNum *u,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)t;
    (void)x;
    (void)u;
    (void)p;
    (void)param;
    (void)userparam;

    // zero matrix placeholder
    for (int i = 0; i < auv::NX * auv::NX; ++i) {
        out[i] = 0.0;
    }
}

void dfdxtrans(typeRNum *out,
               ctypeRNum t,
               ctypeRNum *x,
               ctypeRNum *u,
               ctypeRNum *p,
               const typeGRAMPCparam *param,
               typeUSERPARAM *userparam)
{
    (void)t;
    (void)x;
    (void)u;
    (void)p;
    (void)param;
    (void)userparam;

    // zero matrix placeholder
    for (int i = 0; i < auv::NX * auv::NX; ++i) {
        out[i] = 0.0;
    }
}

void dfdt(typeRNum *out,
          ctypeRNum t,
          ctypeRNum *x,
          ctypeRNum *u,
          ctypeRNum *p,
          const typeGRAMPCparam *param,
          typeUSERPARAM *userparam)
{
    (void)t;
    (void)x;
    (void)u;
    (void)p;
    (void)param;
    (void)userparam;

    // autonomous-system placeholder
    for (int i = 0; i < auv::NX; ++i) {
        out[i] = 0.0;
    }
}

void dHdxdt(typeRNum *out,
            ctypeRNum t,
            ctypeRNum *x,
            ctypeRNum *lambda,
            ctypeRNum *u,
            ctypeRNum *p,
            const typeGRAMPCparam *param,
            typeUSERPARAM *userparam)
{
    (void)t;
    (void)x;
    (void)lambda;
    (void)u;
    (void)p;
    (void)param;
    (void)userparam;

    // autonomous-system placeholder
    for (int i = 0; i < auv::NX; ++i) {
        out[i] = 0.0;
    }
}


} // extern "C"

// ======================================================
// MpcOptim methods
// ======================================================

MpcOptim::MpcOptim() :
    MpcOptim(auv::default_plant_params(),
             auv::default_mpc_weight_params(),
             auv::default_mpc_bounds_params(),
             auv::default_mpc_config())
{
}

MpcOptim::MpcOptim(const auv::AuvPlantParams &plant,
                   const auv::MpcWeightParams &weights,
                   const auv::MpcBoundsParams &bounds,
                   const auv::MpcConfig &config) :
    _plant(plant),
    _weights(weights),
    _bounds(bounds),
    _config(config)
{
    _debug.last_u_nom.zero();
}

int MpcOptim::compute_nhor(double dt, double Thor) const
{
    if (dt <= 0.0) {
        dt = 0.1;
    }

    if (Thor <= 0.0) {
        Thor = 10.0 * dt;
    }

    int nhor = static_cast<int>(floor(Thor / dt + 0.5)) + 1;

    if (nhor < 2) {
        nhor = 2;
    }

    if (nhor > auv::NHOR_MAX) {
        nhor = auv::NHOR_MAX;
    }

    return nhor;
}

void MpcOptim::sync_userparam_from_common()
{
    _userparam.m_u = _plant.m_u;
    _userparam.m_v = _plant.m_v;
    _userparam.Iz  = _plant.Iz;

    _userparam.Du  = _plant.Du;
    _userparam.Dv  = _plant.Dv;
    _userparam.Dr  = _plant.Dr;

    for (int i = 0; i < auv::NX; ++i) {
        _userparam.Q[i] = static_cast<double>(_weights.Q[i]);
        _userparam.P[i] = static_cast<double>(_weights.P[i]);
    }

    for (int i = 0; i < auv::NU; ++i) {
        _userparam.R[i] = static_cast<double>(_weights.R[i]);
    }

    if (_userparam.Nhor_ref < 2) {
        _userparam.Nhor_ref = 2;
    }

    for (int i = 0; i < auv::NX * auv::NHOR_MAX; ++i) {
        _userparam.xref_hor[i] = 0.0;
    }

    for (int i = 0; i < auv::NU * auv::NHOR_MAX; ++i) {
        _userparam.uref_hor[i] = 0.0;
    }

    for (int i = 0; i < auv::NP * auv::NHOR_MAX; ++i) {
        _userparam.dhat_hor[i] = 0.0;
    }
}

void MpcOptim::free_grampc()
{
    if (_grampc != nullptr) {
        grampc_free(&_grampc);
        _grampc = nullptr;
    }

    if (g_up == &_userparam) {
        g_up = nullptr;
    }

    _initialized = false;
    _debug.initialized = false;
}

void MpcOptim::initialize(float dt, float Thor)
{
    free_grampc();
    sync_userparam_from_common();

    double dt_use = (dt > 0.0f) ? static_cast<double>(dt) : 0.02;
    double Thor_use = (Thor > 0.0f) ? static_cast<double>(Thor) : 0.4;

    const int nhor_cfg = (_config.Nhor > 1) ? _config.Nhor : compute_nhor(dt_use, Thor_use);
    const int nhor = (nhor_cfg <= auv::NHOR_MAX) ? nhor_cfg : auv::NHOR_MAX;

    // discrete-time consistency:
    // Thor = (Nhor - 1) * dt
    Thor_use = static_cast<double>(nhor - 1) * dt_use;

    _userparam.Nhor_ref = nhor;
    g_up = &_userparam;

    grampc_init(&_grampc, reinterpret_cast<typeUSERPARAM*>(&_userparam));

    _grampc->param->dt   = dt_use;
    _grampc->param->Thor = Thor_use;
    _grampc->param->t0   = 0.0;

    for (int i = 0; i < auv::NX; ++i) {
        _grampc->param->x0[i]   = 0.0;
        _grampc->param->xdes[i] = 0.0;
    }

    for (int i = 0; i < auv::NU; ++i) {
        _grampc->param->u0[i]   = 0.0;
        _grampc->param->udes[i] = 0.0;
        _grampc->param->umin[i] = static_cast<double>(_bounds.umin[i]);
        _grampc->param->umax[i] = static_cast<double>(_bounds.umax[i]);
    }

    for (int i = 0; i < auv::NP; ++i) {
        _grampc->param->p0[i]   = 0.0;
        _grampc->param->pmin[i] = static_cast<double>(_bounds.pmin[i]);
        _grampc->param->pmax[i] = static_cast<double>(_bounds.pmax[i]);
    }

    grampc_setopt_string(_grampc, "Integrator", "discrete");
    grampc_setopt_string(_grampc, "IntegratorCost", "discrete");
    grampc_setopt_string(_grampc, "LineSearchType", "explicit1");
    grampc_setopt_string(_grampc, "OptimControl", "on");
    grampc_setopt_string(_grampc, "OptimParam", "off");
    grampc_setopt_string(_grampc, "OptimTime", "off");

    grampc_setopt_string(_grampc, "ShiftControl", _config.shift_control ? "on" : "off");
    grampc_setopt_string(_grampc, "TerminalCost", _config.use_terminal_cost ? "on" : "off");
    grampc_setopt_string(_grampc, "IntegralCost", _config.use_integral_cost ? "on" : "off");

    grampc_setopt_int(_grampc, "Nhor", nhor);
    grampc_setopt_int(_grampc, "MaxGradIter", _config.MaxGradIter);
    grampc_setopt_int(_grampc, "MaxMultIter", _config.MaxMultIter);

    _initialized = true;
    _debug.initialized = true;
    _debug.dt = static_cast<float>(dt_use);
    _debug.Thor = static_cast<float>(Thor_use);
    _debug.Nhor = nhor;
    _debug.t0 = 0.0f;
    _debug.last_u_nom.zero();
}

void MpcOptim::set_plant_params(const auv::AuvPlantParams &plant)
{
    _plant = plant;

    if (_initialized) {
        sync_userparam_from_common();
    }
}

void MpcOptim::set_weight_params(const auv::MpcWeightParams &weights)
{
    _weights = weights;

    if (_initialized) {
        sync_userparam_from_common();
    }
}

void MpcOptim::set_bounds_params(const auv::MpcBoundsParams &bounds)
{
    _bounds = bounds;

    if (_initialized && _grampc != nullptr) {
        for (int i = 0; i < auv::NU; ++i) {
            _grampc->param->umin[i] = static_cast<double>(_bounds.umin[i]);
            _grampc->param->umax[i] = static_cast<double>(_bounds.umax[i]);
        }

        for (int i = 0; i < auv::NP; ++i) {
            _grampc->param->pmin[i] = static_cast<double>(_bounds.pmin[i]);
            _grampc->param->pmax[i] = static_cast<double>(_bounds.pmax[i]);
        }
    }
}

void MpcOptim::set_mpc_config(const auv::MpcConfig &config)
{
    _config = config;
}

void MpcOptim::set_reference_horizon(const StateVec &z,
                                     const float *xref_hor,
                                     int nhor,
                                     const ParamVec &dhat)
{
    float dhat_hor[auv::NP * auv::NHOR_MAX]{};

    const int nh = clamp_int(nhor, 2, auv::NHOR_MAX);
    for (int k = 0; k < nh; ++k) {
        for (int i = 0; i < auv::NP; ++i) {
            dhat_hor[k * auv::NP + i] = dhat(i);
        }
    }

    set_reference_horizon_control_preview(z, xref_hor, nullptr, dhat_hor, nh);
}

void MpcOptim::set_reference_horizon_preview(const StateVec &z,
                                             const float *xref_hor,
                                             const float *dhat_hor,
                                             int nhor)
{
    set_reference_horizon_control_preview(z, xref_hor, nullptr, dhat_hor, nhor);
}

void MpcOptim::set_reference_horizon_control_preview(const StateVec &z,
                                                     const float *xref_hor,
                                                     const float *uref_hor,
                                                     const float *dhat_hor,
                                                     int nhor)
{
    if (!_initialized) {
        initialize(0.02f, 0.4f);
    }

    if (_grampc == nullptr) {
        return;
    }

    g_up = &_userparam;

    const int nh = clamp_int(nhor, 2, auv::NHOR_MAX);
    _userparam.Nhor_ref = nh;

    for (int i = 0; i < auv::NX; ++i) {
        _grampc->param->x0[i] = static_cast<double>(z(i));
    }

    if (xref_hor != nullptr) {
        for (int k = 0; k < nh; ++k) {
            for (int i = 0; i < auv::NX; ++i) {
                _userparam.xref_hor[k * auv::NX + i] =
                    static_cast<double>(xref_hor[k * auv::NX + i]);
            }
        }

        for (int k = nh; k < auv::NHOR_MAX; ++k) {
            for (int i = 0; i < auv::NX; ++i) {
                _userparam.xref_hor[k * auv::NX + i] =
                    _userparam.xref_hor[(nh - 1) * auv::NX + i];
            }
        }
    } else {
        for (int k = 0; k < auv::NHOR_MAX; ++k) {
            for (int i = 0; i < auv::NX; ++i) {
                _userparam.xref_hor[k * auv::NX + i] = static_cast<double>(z(i));
            }
        }
    }

    if (uref_hor != nullptr) {
        for (int k = 0; k < nh; ++k) {
            for (int i = 0; i < auv::NU; ++i) {
                _userparam.uref_hor[k * auv::NU + i] =
                    static_cast<double>(uref_hor[k * auv::NU + i]);
            }
        }

        for (int k = nh; k < auv::NHOR_MAX; ++k) {
            for (int i = 0; i < auv::NU; ++i) {
                _userparam.uref_hor[k * auv::NU + i] =
                    _userparam.uref_hor[(nh - 1) * auv::NU + i];
            }
        }
    } else {
        for (int k = 0; k < auv::NHOR_MAX; ++k) {
            for (int i = 0; i < auv::NU; ++i) {
                _userparam.uref_hor[k * auv::NU + i] = 0.0;
            }
        }
    }

    if (dhat_hor != nullptr) {
        for (int k = 0; k < nh; ++k) {
            for (int i = 0; i < auv::NP; ++i) {
                _userparam.dhat_hor[k * auv::NP + i] =
                    static_cast<double>(dhat_hor[k * auv::NP + i]);
            }
        }

        for (int k = nh; k < auv::NHOR_MAX; ++k) {
            for (int i = 0; i < auv::NP; ++i) {
                _userparam.dhat_hor[k * auv::NP + i] =
                    _userparam.dhat_hor[(nh - 1) * auv::NP + i];
            }
        }
    } else {
        for (int k = 0; k < auv::NHOR_MAX; ++k) {
            for (int i = 0; i < auv::NP; ++i) {
                _userparam.dhat_hor[k * auv::NP + i] = 0.0;
            }
        }
    }

    for (int i = 0; i < auv::NX; ++i) {
        _grampc->param->xdes[i] =
            _userparam.xref_hor[(nh - 1) * auv::NX + i];
    }

    for (int i = 0; i < auv::NU; ++i) {
        _grampc->param->udes[i] =
            _userparam.uref_hor[i];
    }

    // 兼容仍然依赖 param->p0 的库内部逻辑：
    // 这里同步写入第 0 个 preview stage。
    for (int i = 0; i < auv::NP; ++i) {
        _grampc->param->p0[i] = _userparam.dhat_hor[i];
    }

    _grampc->param->t0 = static_cast<double>(_debug.t0);
}

MpcOptim::InputVec MpcOptim::solve()
{
    InputVec out{};
    out.zero();

    if (!_initialized || _grampc == nullptr) {
        return out;
    }

    g_up = &_userparam;
    grampc_run(_grampc);

    for (int i = 0; i < auv::NU; ++i) {
        out(i) = static_cast<float>(_grampc->sol->unext[i]);
    }

    _debug.last_u_nom = out;
    _debug.t0 += _debug.dt;

    return out;
}

void MpcOptim::reset()
{
    free_grampc();
    sync_userparam_from_common();
    _debug.last_u_nom.zero();
    _debug.t0 = 0.0f;
}

