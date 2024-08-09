// Gmsh project created on Tue Nov 14 23:47:32 2017
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {2.2, 0, 0, 1.0};
//+
Point(3) = {2.2, 0.41, 0, 1.0};
//+
Point(4) = {0, 0.41, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {0.2, 0.2, 0, 0.05, 0, 2*Pi};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Line Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};

Physical Point(1) = {1, 2, 3, 4};
Physical Line(1) = {1,3,4,5};
Physical Line(2) = {2};
Physical Surface(0) = {1};
