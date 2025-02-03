// Gmsh project created on Wed Oct 16 17:14:54 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, 0, 1, 1.0};
//+
Point(6) = {1, 0, 1, 1.0};
//+
Point(7) = {1, 1, 1, 1.0};
//+
Point(8) = {0, 1, 1, 1.0};
//+
Point(9) = {0.5, 0.5, 1, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Line(9) = {1, 5};
//+
Line(10) = {2, 6};
//+
Line(11) = {3, 7};
//+
Line(12) = {4, 8};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, 10, -5, -9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {3, 12, -7, -11};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 9, -8, -12};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {2, 11, -6, -10};
//+
Plane Surface(6) = {6};

//+
Surface Loop(1) = {1, 3, 6, 4, 5, 2};
//+
Volume(1) = {1};

Physical Point("ptdispxyz", 14) = {4, 1, 2, 3};
//+
Transfinite Curve {4, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12} = 21 Using Progression 1;
//+
Transfinite Surface {1, 2, 3, 4, 5, 6};
Transfinite Volume {:};

Recombine Surface {:};
Recombine Volume {1};

Point{9} In Surface {2};
//Point{9} In Volume {1};
Physical Volume("dom") = {1};
Physical Point("ptforcez", 13) = {9};


// //Mesh the model first to remove duplicate nodes
//Mesh 3;

//Merge "cube.msh";

//// Try optimizing the mesh to avoid duplicate nodes
//Mesh.Optimize = 1;

Mesh 3;

Coherence Mesh;

