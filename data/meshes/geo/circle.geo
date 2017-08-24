// Gmsh project created on Tue Jul 25 13:21:55 2017
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
Circle(1) = {0.5, 0.5, -0, 0.15, 0, 2*Pi};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 1};
//+
Line Loop(1) = {5, 2, 3, 4};
//+
Line Loop(2) = {1};
//+
Plane Surface(1) = {1, 2};

Physical Point(1) = {1, 2, 3, 4};
Physical Line(1) = {1, 2, 4, 5};
Physical Line(2) = {3};
Physical Surface(0) = {1};