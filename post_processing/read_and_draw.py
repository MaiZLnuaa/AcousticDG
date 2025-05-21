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

mesh_file = "../highorder_results/result5000.vtu"
mesh = pv.read(mesh_file)

plotter = pv.Plotter(window_size=[800, 600],title="Pressure [Pa] along X axis")

plotter.add_mesh(mesh, scalars="Pressure [Pa]", show_scalar_bar=False)

plotter.add_scalar_bar( 
                       vertical=True,           # 竖直放置
                       position_x=0.85,         # 横向位置，0~1
                       position_y=0.15,          # 纵向位置，0~1
                       height=0.7,              # 色棒高度比例
                       width=0.1)              # 色棒宽度比例

plotter.view_xy()
plotter.show()

data_name = "Pressure [Pa]"

xmin, xmax = mesh.bounds[0], mesh.bounds[1]
n_samples = 1000

interp = mesh.sample_over_line(pointa=(0, -100, 0), pointb=(100, 0, 0), resolution=n_samples - 1)

x_samples = interp.points[:, 0]
r_samples = np.sqrt(interp.points[:, 0]**2 + (interp.points[:, 1]+100)**2)
y_values = interp[data_name]

print("x_samples.shape:", x_samples.shape)       
print("y_values.shape:", y_values.shape)        
print("Any NaNs in y_values?", np.isnan(y_values).any())

valid = ~np.isnan(y_values)

# 读取 txt 数据
txt_data = np.loadtxt("2d_t100_ma0.5_gassianwall_analytical_data.txt", delimiter=",", )
txt_data = txt_data[txt_data[:, 0].argsort()]
x_txt = txt_data[:, 0]
y_txt = txt_data[:, 1]


plt.plot()
plt.plot(r_samples[valid], y_values[valid], label="Numerical(LEE)", color="blue")
plt.plot(x_txt, y_txt, label="Analytical", color="red", linestyle="dashed")
# plt.xlim([xmin, xmax])
plt.xlim([0, np.max(r_samples)])
plt.ylim([-0.1,0.2])
plt.xlabel("X")
plt.ylabel(data_name)
plt.title(f"{data_name} along X axis")
plt.legend()
plt.grid(True)
plt.show()

