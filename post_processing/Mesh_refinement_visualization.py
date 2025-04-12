#
#  @File          Mesh_refinement_visualization.py
#
#  @Author        MaiZLnuaa <mai-zl@nuaa.edu.cn>
#  @Date          Mon Apr 14 2025 10:48:49
#
#  @Description   python post_processing.py / python post_processing.py --method rbf / python post_processing.py --input result0.vtu --method rbf


import argparse
import meshio
import vtk
import numpy as np
from scipy.interpolate import griddata, RBFInterpolator

# 命令行参数解析
parser = argparse.ArgumentParser(description="Interpolate .vtu data to a denser grid.")
parser.add_argument("--method", choices=["linear", "cubic", "nearest", "rbf"], default="nearest", help="Interpolation method to use")
parser.add_argument("--input", default="result0.vtu", help="Input VTU file")
args = parser.parse_args()

mesh_file = args.input
interp_method = args.method

# Step 1: 读取 .vtu 文件并转为未压缩格式
mesh = meshio.read(mesh_file)
mesh.write("fixed_file.vtu")

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("fixed_file.vtu")
reader.Update()
data = reader.GetOutput()

# Step 2: 判断维度 & 获取原始点和数据
num_points = data.GetNumberOfPoints()
points = np.array([data.GetPoint(i) for i in range(num_points)])

point_data = data.GetPointData().GetArray("Pressure [Pa]")
if point_data is None:
    raise ValueError("❌ 找不到名为 'Pressure [Pa]' 的 PointData！")
values = np.array([point_data.GetValue(i) for i in range(num_points)])

is_2D = np.allclose(points[:, 2], points[0, 2])
print("is 2D" if is_2D else "is 3D")

# Step 3: 生成稠密网格
if is_2D:
    xmin, ymin = points[:, 0].min(), points[:, 1].min()
    xmax, ymax = points[:, 0].max(), points[:, 1].max()
    grid_x, grid_y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    interp_points = np.stack([grid_x.ravel(), grid_y.ravel()], axis=-1)
else:
    xmin, ymin, zmin = points.min(axis=0)
    xmax, ymax, zmax = points.max(axis=0)
    grid_x, grid_y, grid_z = np.mgrid[xmin:xmax:200j, ymin:ymax:200j, zmin:zmax:200j]
    interp_points = np.stack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()], axis=-1)

# Step 4: 插值方法选择
if interp_method == "rbf":
    print("🔵 使用 RBF 插值")
    rbf_interp = RBFInterpolator(points if not is_2D else points[:, :2], values, kernel="thin_plate_spline")
    interp_values = rbf_interp(interp_points)
else:
    method = "cubic" if is_2D else "linear"  
    print("🟡 使用 SciPy griddata 插值, method =",method)
    interp_values = griddata(points if not is_2D else points[:, :2], values, interp_points, method=method) #nearest

interp_values = interp_values.reshape(grid_x.shape)

# Step 5: 创建新的 VTK 网格
output_grid = vtk.vtkUnstructuredGrid()
new_points = vtk.vtkPoints()
new_cells = vtk.vtkCellArray()

# 添加点
if is_2D:
    for i in range(grid_x.shape[0]):
        for j in range(grid_x.shape[1]):
            new_points.InsertNextPoint(grid_x[i, j], grid_y[i, j], 0)
else:
    for i in range(grid_x.shape[0]):
        for j in range(grid_x.shape[1]):
            for k in range(grid_x.shape[2]):
                new_points.InsertNextPoint(grid_x[i, j, k], grid_y[i, j, k], grid_z[i, j, k])

# 添加单元（四边形或六面体）
if is_2D:
    for i in range(grid_x.shape[0] - 1):
        for j in range(grid_x.shape[1] - 1):
            cell = vtk.vtkQuad()
            id0 = i * grid_x.shape[1] + j
            cell.GetPointIds().SetId(0, id0)
            cell.GetPointIds().SetId(1, id0 + 1)
            cell.GetPointIds().SetId(2, id0 + 1 + grid_x.shape[1])
            cell.GetPointIds().SetId(3, id0 + grid_x.shape[1])
            new_cells.InsertNextCell(cell)
else:
    nx, ny, nz = grid_x.shape
    for i in range(nx - 1):
        for j in range(ny - 1):
            for k in range(nz - 1):
                cell = vtk.vtkHexahedron()
                base = i * ny * nz + j * nz + k
                dx = ny * nz
                dy = nz
                cell.GetPointIds().SetId(0, base)
                cell.GetPointIds().SetId(1, base + 1)
                cell.GetPointIds().SetId(2, base + dy + 1)
                cell.GetPointIds().SetId(3, base + dy)
                cell.GetPointIds().SetId(4, base + dx)
                cell.GetPointIds().SetId(5, base + dx + 1)
                cell.GetPointIds().SetId(6, base + dx + dy + 1)
                cell.GetPointIds().SetId(7, base + dx + dy)
                new_cells.InsertNextCell(cell)

# 添加插值结果
pressure_array = vtk.vtkDoubleArray()
pressure_array.SetName("Pressure [Pa]")
pressure_array.SetNumberOfComponents(1)
for v in interp_values.flat:
    pressure_array.InsertNextValue(v)

# 写入新网格
output_grid.SetPoints(new_points)
output_grid.SetCells(9 if is_2D else 12, new_cells)
output_grid.GetPointData().AddArray(pressure_array)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("new_" + mesh_file)
writer.SetInputData(output_grid)
writer.Write()

print(f"✅ 插值完成，输出文件为 new_{mesh_file}")
