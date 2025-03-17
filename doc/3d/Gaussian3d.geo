lc = 5.0;

// 定义立方体的 8 个顶点
Point(1) = {-50, -50, -50, lc};
Point(2) = { 50, -50, -50, lc};
Point(3) = { 50,  50, -50, lc};
Point(4) = {-50,  50, -50, lc};
Point(5) = {-50, -50,  50, lc};
Point(6) = { 50, -50,  50, lc};
Point(7) = { 50,  50,  50, lc};
Point(8) = {-50,  50,  50, lc};

// 定义立方体的 6 个面
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(13) = {1, 2, 3, 4};
Line Loop(14) = {5, 6, 7, 8};
Line Loop(15) = {1, 10, -5, -9};
Line Loop(16) = {2, 11, -6, -10};
Line Loop(17) = {3, 12, -7, -11};
Line Loop(18) = {4, 9, -8, -12};

Plane Surface(19) = {13};
Plane Surface(20) = {14};
Plane Surface(21) = {15};
Plane Surface(22) = {16};
Plane Surface(23) = {17};
Plane Surface(24) = {18};

// 形成体
Surface Loop(25) = {19, 20, 21, 22, 23, 24};
Volume(26) = {25};

Mesh.ElementOrder = 4;  // 生成四阶网格 