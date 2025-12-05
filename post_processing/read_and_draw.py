#
#  @File          read_and_draw.py
#
#  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
#  @Date          Tue May 20 2025 15:25:53
#
#  @Description    
#

import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

this_dir = os.path.dirname(os.path.abspath(__file__))
mesh_file = os.path.join(this_dir, "new_result700.vtu")
# mesh_file = os.path.join(this_dir,  "../highorder_results", "result6000.vtu")
mesh = pv.read(mesh_file)


# # # print(mesh.cell_data)
# # 定义切片平面：Z = 0.0
# z0 = 0.0
# normal = [0, 0, 1]  # Z 方向法向量
# origin = [0, 0, z0]

# # # 创建切片
# slice_mesh = mesh.slice(origin=origin, normal=normal)

# plotter = pv.Plotter(window_size=[800, 600], title=f"Pressure [Pa] on Z={z0} plane")

# plotter.add_mesh(slice_mesh, scalars="Pressure [Pa]", show_scalar_bar=False)



plotter = pv.Plotter(window_size=[800, 600],title="Pressure [Pa] along X axis")

plotter.add_mesh(mesh, scalars="Pressure [Pa]", show_scalar_bar=False)

# plotter.add_scalar_bar( 
#                        vertical=True,           # 竖直放置
#                        position_x=0.85,         # 横向位置，0~1
#                        position_y=0.15,          # 纵向位置，0~1
#                        height=0.7,              # 色棒高度比例
#                        width=0.1)              # 色棒宽度比例

plotter.view_xy()
plotter.show()
plotter.screenshot(os.path.join(this_dir, "pressure_along_x.png"))


data_name = "Pressure [Pa]"

# mesh = mesh.cell_data_to_point_data()  # 将单元数据转换为点数据
# print(mesh.point_data)

# xmin, xmax = mesh.bounds[0], mesh.bounds[1]
n_samples = 10000

interp = mesh.sample_over_line(pointa=(-120, 0, 0), pointb=(120, 0, 0), resolution=n_samples - 1)

x_samples = interp.points[:, 0]
# r_samples = np.sqrt(interp.points[:, 0]**2 + (interp.points[:, 1]+100)**2)
y_values = interp[data_name]

# print("x_samples.shape:", x_samples.shape)       
# print("y_values.shape:", y_values.shape)        
# print("Any NaNs in y_values?", np.isnan(y_values).any())

# valid = ~np.isnan(y_values)

# 读取 txt 数据
# txt_file = os.path.join(this_dir, "2dp4jxj.txt")
# txt_file = os.path.join(this_dir, "2d_t60_ma0.5_gassianwall_analytical_data.txt")
txt_file = os.path.join(this_dir, "2d_t120_ma0.5_gaussian_analytical_data.txt")
# txt_data = np.loadtxt("2d_t120_ma0.5_gaussian_analytical_data.txt", delimiter="," )
# txt_data = np.loadtxt(txt_file)
txt_data = np.loadtxt(txt_file, delimiter=",")
txt_data = txt_data[txt_data[:, 0].argsort()]
x_txt = txt_data[:, 0]
y_txt = txt_data[:, 1]

# output_data = np.column_stack((x_samples, y_values))
# output_txt_file = os.path.join(this_dir, "2dp4jxj.txt")
# np.savetxt(output_txt_file, output_data, fmt="%.8e", delimiter=" ")

plt.plot()
plt.plot(x_samples, y_values, label="Numerical(LEE)", color="blue")
plt.plot(x_txt, y_txt, label="Analytical", color="red", linestyle="dashed")
# plt.xlim([xmin, xmax])
# plt.xlim([-50, np.max(x_samples)])
# plt.ylim([-0.4,0.3])
plt.xlabel("X")
plt.ylabel(data_name)
plt.title(f"{data_name} along X axis")
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(this_dir, "pressure_along_x_data.png"), dpi=300)
plt.show()


# 拟合解析解到插值点：线性插值
from scipy.interpolate import interp1d

# 数值解
x_num = x_samples
u_num = y_values

# 解析解插值到数值解点上
interp_ana = interp1d(x_txt, y_txt, kind="linear", fill_value="extrapolate")
u_exact = interp_ana(x_num)

# L2 误差计算（均匀采样）
dx = (x_num[-1] - x_num[0]) / (len(x_num) - 1)
diff = u_num - u_exact
l2_error = np.sqrt(np.sum(diff**2) * dx)

print(f"L2 error norm: {l2_error:.6e}")


# data_file = os.path.join(this_dir, "..", "results", "observer_1_fft.txt")

# # 读取数据（自动处理列名）
# df = pd.read_csv(data_file, delim_whitespace=True)
# # pd.read_csv(data_file, sep=',')  # 逗号分隔
# # pd.read_csv(data_file, sep=';')  # 分号分隔

# # 画 SPL（Sound Pressure Level） vs Frequency
# plt.figure(figsize=(8, 5))
# plt.plot(df["Frequency"], df["SPL"], color="blue", label="SPL [dB]")
# plt.xlabel("Frequency [Hz or normalized]")
# plt.ylabel("Sound Pressure Level [dB]")
# plt.title("Frequency Spectrum")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig(os.path.join(this_dir, "frequency_spectrum.png"), dpi=300)
# plt.show()
