#
#  @File          3d_section_processing.py
#
#  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
#  @Date          Mon Mar 17 2025 15:28:14
#
#  @Description   Slice the 3D results and interpolate them with a denser grid
#

import meshio
import vtk
import numpy as np
from scipy.interpolate import griddata

mesh_file = "result0.vtu"

# 读取原始 .vtu 文件
mesh = meshio.read(mesh_file)

# 重新保存为未压缩的标准 .vtu
mesh.write("fixed_file.vtu")

# 读取 .vtu 文件
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("fixed_file.vtu")
reader.Update()
data = reader.GetOutput()

# 获取单元的中心坐标
num_cells = data.GetNumberOfCells()
cell_centers = np.zeros((num_cells, 3))

for i in range(num_cells):
    cell = data.GetCell(i)
    cell_points = [data.GetPoint(cell.GetPointId(j)) for j in range(cell.GetNumberOfPoints())]
    cell_centers[i] = np.mean(cell_points, axis=0)  # 计算单元中心

# 选择切面 (例如 z = 0.5)
slice_z = 0.0
tolerance = 0.1  # 允许的误差范围
mask = np.abs(cell_centers[:, 2] - slice_z) < tolerance  # 选择 z ≈ 0.0 的单元

# 只保留切面上的单元
points = cell_centers[mask][:, :2]  # 仅用 x, y 进行插值
values = np.array([data.GetCellData().GetArray("Pressure [Pa]").GetValue(i) for i in range(num_cells)])[mask]

# 生成更密集的网格
xmin, ymin = points[:, 0].min(), points[:, 1].min()
xmax, ymax = points[:, 0].max(), points[:, 1].max()
grid_x, grid_y = np.mgrid[xmin:xmax:300j, ymin:ymax:300j]

# 二维插值
grid_z = griddata(points, values, (grid_x, grid_y), method="cubic")

# 创建 VTK 网格
output_grid = vtk.vtkUnstructuredGrid()
new_points = vtk.vtkPoints()
new_cells = vtk.vtkCellArray()

# 添加点
for i in range(grid_x.shape[0]):
    for j in range(grid_x.shape[1]):
        new_points.InsertNextPoint(grid_x[i, j], grid_y[i, j], slice_z)  # 仍然保持原来的 z

# 创建单元（四边形单元）
for i in range(grid_x.shape[0] - 1):
    for j in range(grid_x.shape[1] - 1):
        cell = vtk.vtkQuad()
        cell.GetPointIds().SetId(0, i * grid_x.shape[1] + j)
        cell.GetPointIds().SetId(1, i * grid_x.shape[1] + (j + 1))
        cell.GetPointIds().SetId(2, (i + 1) * grid_x.shape[1] + (j + 1))
        cell.GetPointIds().SetId(3, (i + 1) * grid_x.shape[1] + j)
        new_cells.InsertNextCell(cell)

# 处理插值后的数据
pressure_array = vtk.vtkDoubleArray()
pressure_array.SetName("Pressure [Pa]")

for i in range(grid_z.size):
    pressure_array.InsertNextValue(grid_z.flat[i])

# 组合数据
output_grid.SetPoints(new_points)
output_grid.SetCells(9, new_cells)  # 四边形单元
output_grid.GetPointData().AddArray(pressure_array)

# 写入 .vtu 文件
output_file = "slice_" + mesh_file
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(output_file)
writer.SetInputData(output_grid)
writer.Write()
print(f"✅ 切面数据写入完成，生成 {output_file} ！")
