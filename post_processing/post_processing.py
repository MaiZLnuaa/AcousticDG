#
#  @File          post_processing.py
#
#  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
#  @Date          Mon Mar 14 2025 15:28:42
#
#  @Description   Interpolate the results using a denser grid
#


import meshio
import vtk
import numpy as np
from scipy.interpolate import griddata

mesh_file = "result1200_2d_23246_p2.vtu"

# 读取原始 .vtu 文件
mesh = meshio.read(mesh_file)

# # 打印基本信息
# print(mesh)

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

# 判断数据是 2D 还是 3D
z_values = cell_centers[:, 2]
is_2D = np.allclose(z_values, z_values[0])  # 如果所有 z 值都相等，说明是 2D

# 选择适当的坐标进行插值
if is_2D:
    points = cell_centers[:, :2]  # 仅用 x, y
else:
    points = cell_centers  # 直接使用 (x, y, z)

# 获取某个标量场数据
cell_data = data.GetCellData().GetArray("Pressure [Pa]")
values = np.array([cell_data.GetValue(i) for i in range(cell_data.GetNumberOfTuples())])

# 生成更密集的网格
if is_2D:
    xmin, ymin = points[:, 0].min(), points[:, 1].min()
    xmax, ymax = points[:, 0].max(), points[:, 1].max()
    grid_x, grid_y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    
    # 二维插值
    grid_z = griddata(points, values, (grid_x, grid_y), method="cubic")

else:
    xmin, ymin, zmin = points[:, 0].min(), points[:, 1].min(), points[:, 2].min()
    xmax, ymax, zmax = points[:, 0].max(), points[:, 1].max(), points[:, 2].max()
    grid_x, grid_y, grid_z = np.mgrid[xmin:xmax:100j, ymin:ymax:100j, zmin:zmax:100j]

    # 三维插值
    grid_values = griddata(points, values, (grid_x, grid_y, grid_z), method="nearest") #method="nearest"

# 创建 VTK 网格
output_grid = vtk.vtkUnstructuredGrid()
new_points = vtk.vtkPoints()
new_cells = vtk.vtkCellArray()

# 添加点
if is_2D:
    for i in range(grid_x.shape[0]):
        for j in range(grid_x.shape[1]):
            new_points.InsertNextPoint(grid_x[i, j], grid_y[i, j], 0)  # 2D：z = 0
else:
    for i in range(grid_x.shape[0]):
        for j in range(grid_x.shape[1]):
            for k in range(grid_x.shape[2]):
                new_points.InsertNextPoint(grid_x[i, j, k], grid_y[i, j, k], grid_z[i, j, k])  # 3D

# 创建单元（四边形或六面体）
if is_2D:
    for i in range(grid_x.shape[0] - 1):
        for j in range(grid_x.shape[1] - 1):
            cell = vtk.vtkQuad()
            cell.GetPointIds().SetId(0, i * grid_x.shape[1] + j)
            cell.GetPointIds().SetId(1, i * grid_x.shape[1] + (j + 1))
            cell.GetPointIds().SetId(2, (i + 1) * grid_x.shape[1] + (j + 1))
            cell.GetPointIds().SetId(3, (i + 1) * grid_x.shape[1] + j)
            new_cells.InsertNextCell(cell)
else:
    for i in range(grid_x.shape[0] - 1):
        for j in range(grid_x.shape[1] - 1):
            for k in range(grid_x.shape[2] - 1):
                cell = vtk.vtkHexahedron()
                base_index = i * grid_x.shape[1] * grid_x.shape[2] + j * grid_x.shape[2] + k
                cell.GetPointIds().SetId(0, base_index)
                cell.GetPointIds().SetId(1, base_index + 1)
                cell.GetPointIds().SetId(2, base_index + grid_x.shape[2] + 1)
                cell.GetPointIds().SetId(3, base_index + grid_x.shape[2])
                cell.GetPointIds().SetId(4, base_index + grid_x.shape[1] * grid_x.shape[2])
                cell.GetPointIds().SetId(5, base_index + grid_x.shape[1] * grid_x.shape[2] + 1)
                cell.GetPointIds().SetId(6, base_index + grid_x.shape[1] * grid_x.shape[2] + grid_x.shape[2] + 1)
                cell.GetPointIds().SetId(7, base_index + grid_x.shape[1] * grid_x.shape[2] + grid_x.shape[2])
                new_cells.InsertNextCell(cell)

# 处理插值后的数据
pressure_array = vtk.vtkDoubleArray()
pressure_array.SetName("Pressure [Pa]")

if is_2D:
    for i in range(grid_z.size):
        pressure_array.InsertNextValue(grid_z.flat[i])
else:
    for i in range(grid_values.size):
        pressure_array.InsertNextValue(grid_values.flat[i])

# 组合数据
output_grid.SetPoints(new_points)
if is_2D:
    output_grid.SetCells(9, new_cells)  # 四边形单元
else:
    output_grid.SetCells(12, new_cells)  # 六面体单元
output_grid.GetPointData().AddArray(pressure_array)

# 写入 .vtu 文件
output_file = "new_" + mesh_file
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(output_file)
writer.SetInputData(output_grid)
writer.Write()
print(f"✅ 写入完成，生成 {output_file} ！")
