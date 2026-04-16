uint64 timestamp          # 系统启动以来的时间 (microseconds)

float32 sim_time          # 当前仿真时间 (s)
float32 dt                # 当前实际使用的时间步长 (s)

# 核心状态变量 [X, Y, psi, u, v, r]
float32[6] x              # Actual state (真实状态)
float32[6] z              # Nominal state (名义状态)
float32[6] xref           # Reference state (当前上层任务参考点)
float32[6] error_x_z      # Tube误差 (x - z)

# 扰动项 [d_u, d_v, d_r]
float32[3] d_true         # 真实外扰 (来自 WaveDisturbance)
float32[3] d_hat          # 学习外扰 (来自 HopfObserver)

# 控制量 [tau_u, tau_v, tau_r]
# 控制量与执行器调试 [tau_u, tau_v, tau_r]
float32[3] tau_nom        # 第一层 MPC 计算的名义控制量
float32[3] tau_fb         # 闭环反馈控制量 (原代码预留)
float32[3] tau_ancillary  # 第二层 MPC 计算的辅助控制量
float32[3] tau_raw        # 限幅前的原始指令 (tau_cmd)
float32[3] tau_applied    # 最终施加到系统的控制量 (含限幅)
float32[3] tau_clip       # 控制量被物理限幅截断的部分

# Hopf 观测器详细调试 [d_raw, d_pred_before, f_error]
float32[3] hopf_surge
float32[3] hopf_sway
float32[3] hopf_yaw

# --- 新增：算法性能评估指标 (用于 PlotJuggler 监控和日志分析) ---
float32 disturbance_residual_norm    # 扰动学习残差范数 (resid)
float32 preview_update_norm          # 预测域更新偏差范数 (prvupd)
float32 preview_span_norm            # 预测域内波动跨度 (prvspn)
float32 tube_tracking_error_norm     # Tube 跟踪误差总范数 (tubee)
float32 tube_value_surrogate         # 当前 Lyapunov/代价函数值 (tubev)
float32 tube_value_delta             # 代价函数变化率 (tubedv)

# --- 新增：约束违反统计 (正值表示越界) ---
float32 central_state_violation      # 名义层状态违反范数 (zviol)
float32 anc_state_violation          # 真实层状态违反范数 (xviol)
float32 central_input_violation      # 名义层输入越界范数 (cvu)
float32 anc_input_violation          # 辅助层输入越界范数 (avu)

# --- 新增：逻辑标志位 (Bool 在 msg 中通常映射为 uint8) ---
bool hopf_preview_ok                 # Hopf 预测域外推是否成功
bool central_export_ok               # 第一层时域解导出是否成功
bool central_bridge_ok               # GRAMPC 实例 1 桥接是否匹配
bool anc_bridge_ok                   # GRAMPC 实例 2 桥接是否匹配

uint8 xref_nhor                      # 当前预测域实际步数
