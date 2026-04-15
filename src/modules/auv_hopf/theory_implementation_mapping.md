# Hopf + Tube-MPC 理论—实现映射

这份文件把论文里的关键理论对象，和当前代码中的真实对象一一对齐，供整改验收与答辩使用。

## 两层控制结构

- 第一层 central path / nominal OCP  
  代码对象：`AuvHopf::_central_mpc`，输入参考为 `_mission_hor`，输出为 `_v_central`，并导出 `_zref_hor / _vref_hor`。
- 第二层 Hopf-augmented ancillary NMPC  
  代码对象：`AuvHopf::_anc_mpc`，状态参考 `_zref_hor`，控制参考 `_vref_hor`，扰动预览 `_dhat_hor`。
- 执行层 physical saturation  
  代码对象：`AuvHopf::_actuator_sat`。

## 扰动与 preview

- 真扰动 `d_k`  
  代码对象：`AuvHopf::_d_true`
- Hopf 学习估计 `d_hat[k|k]`  
  代码对象：`AuvHopf::_d_hat`
- 多步预览 `d_hat[k+i|k]`  
  代码对象：`AuvHopf::_dhat_hor`
  语义：由当前已学得 Hopf 内模状态做 short-horizon signal extrapolation 得到的外生 preview，
  不是基于未来海浪真参数与未来状态耦合关系得到的真扰动预报。

## practical ISS / tube 相关在线量

- 残余扰动项 `\bar d_r` 的在线代理  
  代码对象：`AuvHopf::_disturbance_residual_norm = ||d_true - d_hat||`
- preview-update 项 `\Delta d` 的在线代理  
  代码对象：`AuvHopf::_preview_update_norm`
- 相对 central path 的 tracking error `||x_k - z_k^*||`  
  代码对象：`AuvHopf::_tube_tracking_error_norm`
- ancillary value function 的在线二次型代理  
  代码对象：`AuvHopf::_tube_value_surrogate`
- 一步 value drift 代理  
  代码对象：`AuvHopf::_tube_value_delta`

## 约束分层

- nominal / central 输入边界  
  代码对象：`_central_mpc.bounds_params()`
- ancillary 输入边界  
  代码对象：`_anc_mpc.bounds_params()`
- central 状态边界  
  代码对象：`AuvHopf::_central_state_bounds`
- ancillary 状态边界  
  代码对象：`AuvHopf::_anc_state_bounds`
- physical actuator limits  
  代码对象：`AuvHopf::_physical_limits` 与 `ActuatorSaturation`

## 在线验证指标

- nominal 状态越界量  
  代码对象：`_central_state_violation_norm`
- ancillary / actual 状态越界量  
  代码对象：`_anc_state_violation_norm`
- central 输入越界量  
  代码对象：`_central_input_violation_norm`
- ancillary 输入越界量  
  代码对象：`_anc_input_violation_norm`
- physical clipping 量  
  代码对象：`_physical_clip_norm`

## GRAMPC 边界条件防护

- 当前桥接 owner 是否匹配当前实例  
  代码对象：`MpcOptim::DebugInfo::bridge_owner_match`
- 最近一次时域导出是否成功  
  代码对象：`MpcOptim::DebugInfo::last_export_ok`
- RODAS 占位函数是否被误触发  
  代码对象：`MpcOptim::DebugInfo::rodas_stub_touch_count`
