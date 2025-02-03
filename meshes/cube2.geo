// Define points for the cube geometry
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {0, 0, 1, 1.0};
Point(6) = {1, 0, 1, 1.0};
Point(7) = {1, 1, 1, 1.0};
Point(8) = {0, 1, 1, 1.0};

// Define an additional internal point that needs to be in a physical group
Point(100) = {0.5, 0.5, 0.5, 1.0};  // This is the internal point

// Define lines for the cube
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

// Define surfaces for the cube
Line Loop(13) = {1, 2, 3, 4};
Plane Surface(14) = {13};
Line Loop(15) = {5, 6, 7, 8};
Plane Surface(16) = {15};
Line Loop(17) = {1, 10, -5, -9};
Plane Surface(18) = {17};
Line Loop(19) = {2, 11, -6, -10};
Plane Surface(20) = {19};
Line Loop(21) = {3, 12, -7, -11};
Plane Surface(22) = {21};
Line Loop(23) = {4, 9, -8, -12};
Plane Surface(24) = {23};

// Transfinite lines and surfaces for uniform mesh
Transfinite Line {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12} = 3;  // 3 divisions per line
Transfinite Surface {14, 16, 18, 20, 22, 24};

// Define the volume for the cube
Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25};
Transfinite Volume {26};

// Recombine for hexahedral mesh
Recombine Volume {26};

// Define the internal point and embed it into the mesh
Point{100} In Volume {26}; // This embeds the point into the volume mesh

// Physical group for the volume
Physical Volume("CubeVolume") = {26};

// Physical group for the internal point
Physical Point("InternalPoint") = {100};

// Optimize the mesh to avoid duplicate nodes
Mesh.Optimize = 1;
