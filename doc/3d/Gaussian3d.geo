// lc = 5.0;
// half = 50;
// // 定义立方体的 8 个顶点
// Point(1) = {-half, -half, -half, lc};
// Point(2) = { half, -half, -half, lc};
// Point(3) = { half,  half, -half, lc};
// Point(4) = {-half,  half, -half, lc};
// Point(5) = {-half, -half,  half, lc};
// Point(6) = { half, -half,  half, lc};
// Point(7) = { half,  half,  half, lc};
// Point(8) = {-half,  half,  half, lc};

// // 定义边
// Line(1) = {1, 2};
// Line(2) = {2, 3};
// Line(3) = {3, 4};
// Line(4) = {4, 1};
// Line(5) = {5, 6};
// Line(6) = {6, 7};
// Line(7) = {7, 8};
// Line(8) = {8, 5};
// Line(9) = {1, 5};
// Line(10) = {2, 6};
// Line(11) = {3, 7};
// Line(12) = {4, 8};

// Line Loop(13) = {1, 2, 3, 4};
// Line Loop(14) = {5, 6, 7, 8};
// Line Loop(15) = {1, 10, -5, -9};
// Line Loop(16) = {2, 11, -6, -10};
// Line Loop(17) = {3, 12, -7, -11};
// Line Loop(18) = {4, 9, -8, -12};

// Plane Surface(19) = {13};
// Plane Surface(20) = {14};
// Plane Surface(21) = {15};
// Plane Surface(22) = {16};
// Plane Surface(23) = {17};
// Plane Surface(24) = {18};

// Physical Surface("Absorbing") = {13, 14, 15, 16, 17, 18};

// // 形成体
// Surface Loop(25) = {19, 20, 21, 22, 23, 24};
// Volume(26) = {25};
// physical Volume("Domain") = {26};


// Mesh.ElementOrder = 3;  // 生成四阶网格 




lc = 5.0;
half = 50;

// 上半点 (z = 0 到 +50)
Point(1) = {-half, -half,  0, lc};
Point(2) = { half, -half,  0, lc};
Point(3) = { half,  half,  0, lc};
Point(4) = {-half,  half,  0, lc};
Point(5) = {-half, -half,  half, lc};
Point(6) = { half, -half,  half, lc};
Point(7) = { half,  half,  half, lc};
Point(8) = {-half,  half,  half, lc};

// 下半点 (z = 0 到 -50)
Point(9)  = {-half, -half, -half, lc};
Point(10) = { half, -half, -half, lc};
Point(11) = { half,  half, -half, lc};
Point(12) = {-half,  half, -half, lc};

// 上半部的面
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
Line(9) = {1,5}; Line(10) = {2,6}; Line(11) = {3,7}; Line(12) = {4,8};

Line Loop(13) = {1,2,3,4};               Plane Surface(19) = {13}; // bottom @ z=0
Line Loop(14) = {5,6,7,8};               Plane Surface(20) = {14}; // top @ z=+50
Line Loop(15) = {1,10,-5,-9};            Plane Surface(21) = {15};
Line Loop(16) = {2,11,-6,-10};           Plane Surface(22) = {16};
Line Loop(17) = {3,12,-7,-11};           Plane Surface(23) = {17};
Line Loop(18) = {4,9,-8,-12};            Plane Surface(24) = {18};

Surface Loop(25) = {19,20,21,22,23,24};  
Volume(26) = {25};                       // 上半体

// 下半部边和面
Line(101) = {9,10}; Line(102) = {10,11}; Line(103) = {11,12}; Line(104) = {12,9};
Line(105) = {9,1}; Line(106) = {10,2}; Line(107) = {11,3}; Line(108) = {12,4};

Line Loop(109) = {101,102,103,104};      Plane Surface(119) = {109}; // bottom @ z=-50
Line Loop(110) = {-1,-2,-3,-4};          Plane Surface(120) = {110}; // top @ z=0 (reuse)
Line Loop(111) = {101,106,-1,-105};      Plane Surface(121) = {111};
Line Loop(112) = {102,107,-2,-106};      Plane Surface(122) = {112};
Line Loop(113) = {103,108,-3,-107};      Plane Surface(123) = {113};
Line Loop(114) = {104,105,-4,-108};      Plane Surface(124) = {114};

Surface Loop(125) = {119,120,121,122,123,124};
Volume(126) = {125};                      // 下半体

// 物理组
Physical Volume("TopHalf") = {26};
Physical Volume("BottomHalf") = {126};
Physical Surface("Absorbing") = {19,20,21,22,23,24,119,120,121,122,123,124};

Mesh.ElementOrder = 3;



