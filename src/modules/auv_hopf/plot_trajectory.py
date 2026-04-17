import re
import matplotlib.pyplot as plt
import numpy as np


def plot_hpp_trajectory(filename):
    t_list, x_list, y_list, psi_list = [], [], [], []
    u_list, v_list, r_list = [], [], []
    tau_u_list, tau_v_list, tau_r_list = [], [], []

    print(f"正在解析 {filename} ...")

    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            match = re.search(r'\{([^}]+)\}', line)
            if not match:
                continue

            data_str = match.group(1).replace('f', '')
            try:
                data = [float(val.strip()) for val in data_str.split(',') if val.strip()]
            except ValueError:
                continue

            # 新格式: {t, X, Y, psi, u, v, r, tau_u, tau_v, tau_r}
            if len(data) == 10:
                t_list.append(data[0])
                x_list.append(data[1])
                y_list.append(data[2])
                psi_list.append(data[3])
                u_list.append(data[4])
                v_list.append(data[5])
                r_list.append(data[6])
                tau_u_list.append(data[7])
                tau_v_list.append(data[8])
                tau_r_list.append(data[9])

    if not t_list:
        print("未能在文件中找到有效的轨迹数据！")
        return

    print(f"解析成功！共提取到 {len(t_list)} 个轨迹点。正在绘图...")

    t_arr = np.array(t_list)
    x_arr = np.array(x_list)
    y_arr = np.array(y_list)
    psi_arr = np.array(psi_list)
    u_arr = np.array(u_list)
    v_arr = np.array(v_list)
    r_arr = np.array(r_list)
    tau_u_arr = np.array(tau_u_list)
    tau_v_arr = np.array(tau_v_list)
    tau_r_arr = np.array(tau_r_list)

    fig = plt.figure(figsize=(15, 10))

    # 1. 2D 轨迹
    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(y_arr, x_arr, label='Reference Trajectory', color='b', linewidth=2)
    ax1.scatter(y_arr[0], x_arr[0], color='g', marker='o', s=80, label='Start')
    ax1.scatter(y_arr[-1], x_arr[-1], color='r', marker='x', s=80, label='End')


    ax1.set_title('AUV 2D Reference Path (NED Frame)')
    ax1.set_xlabel('Y (East) [m]')
    ax1.set_ylabel('X (North) [m]')
    ax1.axis('equal')
    ax1.grid(True, linestyle='--')
    ax1.legend()

    # 2. 速度曲线
    ax2 = plt.subplot(2, 2, 2)
    ax2.plot(t_arr, u_arr, label='Surge u', color='b')
    ax2.plot(t_arr, v_arr, label='Sway v', color='r')
    ax2.plot(t_arr, r_arr, label='Yaw Rate r', color='purple')
    ax2.set_title('Velocities')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('State')
    ax2.grid(True)
    ax2.legend()

    # 3. 偏航角
    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(t_arr, psi_arr, label='Yaw angle psi', color='darkorange')
    ax3.set_title('Yaw Angle')
    ax3.set_xlabel('Time [s]')
    ax3.set_ylabel('rad')
    ax3.grid(True)
    ax3.legend()

    # 4. 标称控制输入
    ax4 = plt.subplot(2, 2, 4)
    ax4.plot(t_arr, tau_u_arr, label='tau_u', color='b')
    ax4.plot(t_arr, tau_v_arr, label='tau_v', color='r')
    ax4.plot(t_arr, tau_r_arr, label='tau_r', color='purple')
    ax4.set_title('Nominal Inputs')
    ax4.set_xlabel('Time [s]')
    ax4.set_ylabel('Input')
    ax4.grid(True)
    ax4.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_hpp_trajectory("trajectory.hpp")