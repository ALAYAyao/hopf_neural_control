uint64 timestamp          # time since system start (microseconds)

float32 sim_time          # 当前仿真时间 (s)
float32 dt                # 当前实际使用的时间步长 (s)

# 状态变量 [X, Y, psi, u, v, r]
float32[6] x              # Actual state (真实状态)
float32[6] z              # Nominal state (名义状态)
float32[6] xref           # Reference state (当前参考点)
float32[6] error_x_z      # Tube误差 (x - z)

# 扰动项 [d_u, d_v, d_r]
float32[3] d_true         # 真实外扰 (来自 WaveDisturbance)
float32[3] d_hat          # 学习外扰 (来自 HopfObserver)

# 控制量 [tau_u, tau_v, tau_r]
float32[3] tau_nom        # MPC 计算的名义控制量
float32[3] tau_applied    # 最终施加的控制量 (含 Tube 补偿)

# Hopf 观测器调试信息 [d_raw, d_pred_before, f_error]
float32[3] hopf_surge
float32[3] hopf_sway
float32[3] hopf_yaw

# 执行器分配/辅助控制调试 [tau_u, tau_v, tau_r]
float32[3] tau_fb
float32[3] tau_ancillary
float32[3] tau_raw

uint8 xref_nhor           # 当前预测域实际步数
