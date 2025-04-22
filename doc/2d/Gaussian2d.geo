lc = 2.0;
half = 10.0;
Point(1) = {-half, -half, 0, lc};
Point(2) = { half, -half, 0, lc};
Point(3) = { half,  half, 0, lc};
Point(4) = {-half,  half, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// 物面边界（Reflecting）：下边界（1）+ 右边界（2）
Physical Curve("Reflecting") = {1, 2};

// 远场边界（Absorbing）：左边界（4）+ 上边界（3）
Physical Curve("Absorbing") = {3, 4};

// 计算域
Physical Surface("Domain") = {1};

Mesh.ElementOrder = 2;  // 生成四阶网格

Mesh 2;
