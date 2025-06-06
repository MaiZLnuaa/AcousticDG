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

this_dir = os.path.dirname(os.path.abspath(__file__))
# mesh_file = os.path.join(this_dir, "result1000.vtu")
mesh_file = os.path.join(this_dir, "..", "highorder_results", "result1200.vtu")
mesh = pv.read(mesh_file)


# # print(mesh.cell_data)
# # 定义切片平面：Z = 0.0
# z0 = 0.0
# normal = [0, 0, 1]  # Z 方向法向量
# origin = [0, 0, z0]

# # # 创建切片
# slice_mesh = mesh.slice(origin=origin, normal=normal)

# plotter = pv.Plotter(window_size=[800, 600], title=f"Pressure [Pa] on Z={z0} plane")

# plotter.add_mesh(slice_mesh, scalars="p_prime", show_scalar_bar=False)



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

data_name = "Pressure [Pa]"

# mesh = mesh.cell_data_to_point_data()  # 将单元数据转换为点数据
# print(mesh.point_data)

# xmin, xmax = mesh.bounds[0], mesh.bounds[1]
n_samples = 2000

interp = mesh.sample_over_line(pointa=(-50, 0, 0), pointb=(50, 0, 0), resolution=n_samples - 1)

x_samples = interp.points[:, 0]
# r_samples = np.sqrt(interp.points[:, 0]**2 + (interp.points[:, 1]+100)**2)
y_values = interp[data_name]

# print("x_samples.shape:", x_samples.shape)       
# print("y_values.shape:", y_values.shape)        
# print("Any NaNs in y_values?", np.isnan(y_values).any())

# valid = ~np.isnan(y_values)

# 读取 txt 数据
txt_file = os.path.join(this_dir, "monopole_t270.txt")
# txt_data = np.loadtxt("2d_t60_ma0.5_gassianwall_analytical_data.txt", delimiter=",", )
txt_data = np.loadtxt(txt_file)
txt_data = txt_data[txt_data[:, 0].argsort()]
x_txt = txt_data[:, 0]
y_txt = txt_data[:, 1]

# output_data = np.column_stack((x_samples, y_values))
# output_txt_file = os.path.join(this_dir, "2dp4jxj.txt")
# np.savetxt(output_txt_file, output_data, fmt="%.8e", delimiter=" ")

plt.plot()
plt.plot(x_samples, y_values, label="Numerical(LEE)", color="blue")
# plt.plot(x_txt, y_txt, label="Analytical", color="red", linestyle="dashed")
# plt.xlim([xmin, xmax])
# plt.xlim([-50, np.max(x_samples)])
# plt.ylim([-0.4,0.3])
plt.xlabel("X")
plt.ylabel(data_name)
plt.title(f"{data_name} along X axis")
plt.legend()
plt.grid(True)
plt.show()

