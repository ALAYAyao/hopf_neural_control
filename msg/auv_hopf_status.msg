uint64 timestamp          # 系统启动以来的时间 (microseconds)

float32 sim_time          # 当前仿真时间 (s)
float32 dt                # 当前实际使用的时间步长 (s)

# 核心状态变量 [X, Y, psi, u, v, r]
float32[6] x              # Actual state (真实状态)
float32[6] z              # Nominal state (优化/名义轨迹)
float32[6] xref           # Reference state (跟踪参考轨迹)
float32[6] error_x_z      # 跟踪残差 (x - z)

# 扰动项 [d_u, d_v, d_r]
float32[3] d_true         # 真实外扰
float32[3] d_hat          # 学习外扰
float32[3] d_res          # 干扰学习残差 (d_true - d_hat)，按 surge/sway/yaw 三通道给出

float32[3] tau_applied          #surge/sway/yaw 三通道给出
# 核心标量指标
float32 tube_tracking_error_norm   # 跟踪残差总范数
