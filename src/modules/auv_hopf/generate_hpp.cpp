#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// ----------------------------------------
// 角度归一化到 [-pi, pi]
// 用 fmodf 比 while 循环更稳更快
// ----------------------------------------
static inline float wrap_pi(float angle)
{
    constexpr float TWO_PI = 2.0f * static_cast<float>(M_PI);
    angle = std::fmod(angle + static_cast<float>(M_PI), TWO_PI);
    if (angle < 0.0f) {
        angle += TWO_PI;
    }
    return angle - static_cast<float>(M_PI);
}

// ----------------------------------------
// 参考轨迹一个采样点
// Data format:
// {t, X, Y, psi, u, v, r, tau_u, tau_v, tau_r}
// ----------------------------------------
struct RefPoint {
    float t;
    float X;
    float Y;
    float psi;
    float u;
    float v;
    float r;
    float tau_u;
    float tau_v;
    float tau_r;
};

int main()
{
    // ========================================
    // 1. 名义 AUV 参数
    // 对应当前 default_plant_params()
    // ========================================
    constexpr float m_u = 25.0f;
    constexpr float m_v = 35.0f;
    constexpr float Iz  = 3.5f;
    constexpr float Du  = 12.0f;
    constexpr float Dv  = 25.0f;
    constexpr float Dr  = 8.0f;

    // ========================================
    // 2. 标称输入 tau_d
    // 这组值决定离线参考轨迹形状
    // ========================================
    constexpr float tau_u = 6.0f;
    constexpr float tau_v = 0.0f;
    constexpr float tau_r = 0.2f;

    // ========================================
    // 3. 初始状态 x_d(0) = [X,Y,psi,u,v,r]
    // ========================================
    float X   = 20.0f;
    float Y   = 20.0f;
    float psi = static_cast<float>(M_PI / 2.0);
    float u   = 0.0f;
    float v   = 0.0f;
    float r   = 0.0f;

    // ========================================
    // 4. 离线仿真参数
    // ========================================
    constexpr float dt    = 0.1f;
    constexpr float T_max = 150.0f;
    constexpr int TRAJ_LENGTH = static_cast<int>(T_max / dt) + 1;

    // ========================================
    // 5. 输出头文件
    // ========================================
    std::ofstream outfile("trajectory.hpp", std::ios::out | std::ios::trunc);
    if (!outfile.is_open()) {
        std::cerr << "打开 trajectory.hpp 失败！" << std::endl;
        return -1;
    }

    outfile << "#pragma once\n\n";
    outfile << "namespace auv {\n\n";
    outfile << "static constexpr int TRAJ_LENGTH = " << TRAJ_LENGTH << ";\n\n";
    outfile << "// Data format: {t, X, Y, psi, u, v, r, tau_u, tau_v, tau_r}\n";
    outfile << "static const float REF_TRAJECTORY[TRAJ_LENGTH][10] = {\n";

    outfile << std::fixed << std::setprecision(6);

    // ========================================
    // 6. 直接离线积分 + 直接写文件
    // 不做大数组缓存，节省内存
    // ========================================
    for (int k = 0; k < TRAJ_LENGTH; ++k) {
        const float t = static_cast<float>(k) * dt;

        // 写当前点
        outfile << "    {"
                << t     << "f, "
                << X     << "f, "
                << Y     << "f, "
                << psi   << "f, "
                << u     << "f, "
                << v     << "f, "
                << r     << "f, "
                << tau_u << "f, "
                << tau_v << "f, "
                << tau_r << "f}";

        if (k < TRAJ_LENGTH - 1) {
            outfile << ",";
        }
        outfile << "\n";

        // 最后一个点只写，不再推进
        if (k == TRAJ_LENGTH - 1) {
            break;
        }

        // ----------------------------------------
        // 运动学
        // ----------------------------------------
        const float X_dot   = std::cos(psi) * u - std::sin(psi) * v;
        const float Y_dot   = std::sin(psi) * u + std::cos(psi) * v;
        const float psi_dot = r;

        // ----------------------------------------
        // 动力学（无扰动）
        // 与当前工程名义模型一致
        // ----------------------------------------
        const float u_dot = (tau_u - Du * u + m_v * v * r) / m_u;
        const float v_dot = (tau_v - Dv * v - m_u * u * r) / m_v;
        const float r_dot = (tau_r - Dr * r + (m_u - m_v) * u * v) / Iz;

        // ----------------------------------------
        // 欧拉步进
        // ----------------------------------------
        X   += dt * X_dot;
        Y   += dt * Y_dot;
        psi  = wrap_pi(psi + dt * psi_dot);
        u   += dt * u_dot;
        v   += dt * v_dot;
        r   += dt * r_dot;
    }

    outfile << "};\n\n";
    outfile << "} // namespace auv\n";
    outfile.close();

    std::cout << "成功生成 trajectory.hpp，共 " << TRAJ_LENGTH << " 个轨迹点。" << std::endl;
    return 0;
}