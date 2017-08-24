// Gmsh project created on Tue Jul 25 13:00:20 2017
h = 10.0;
//+
SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, h};
//+
Point(2) = {3, 0, 0, h};
//+
Point(3) = {3, 1, 0, h};
//+
Point(4) = {0, 1, 0, h};
//+
Circle(1) = {0.5, 0.7, 0, 0.2, 0, 2*Pi};
//+
Circle(2) = {0.3, 0.3, 0, 0.1, 0, 2*Pi};
//+
Circle(3) = {1.1, 0.4, 0, 0.3, 0, 2*Pi};
//+
Circle(4) = {2.5, 0.7, 0, 0.2, 0, 2*Pi};
//+
Circle(5) = {1.5, 0.8, 0, 0.1, 0, 2*Pi};
//+
Circle(6) = {1.7, 0.4, 0, 0.1, 0, 2*Pi};
//+
Circle(7) = {1.9, 0.7, 0, 0.1, 0, 2*Pi};
//+
Circle(8) = {1.5, 0.2, 0, 0.1, 0, 2*Pi};
//+
Circle(9) = {2, 0.3, 0, 0.15, 0, 2*Pi};
//+
Circle(10) = {2.5, 0.2, 0, 0.15, 0, 2*Pi};
//+
Line(11) = {1, 2};
//+
Line(12) = {2, 3};
//+
Line(13) = {3, 4};
//+
Line(14) = {4, 1};
//+
Line Loop(1) = {14, 11, 12, 13};
//+
Line Loop(2) = {1};
//+
Line Loop(3) = {2};
//+
Line Loop(4) = {3};
//+
Line Loop(5) = {4};
//+
Line Loop(6) = {5};
//+
Line Loop(7) = {6};
//+
Line Loop(8) = {7};
//+
Line Loop(9) = {8};
//+
Line Loop(10) = {9};
//+
Line Loop(11) = {10};
//+
Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

Physical Point(1) = {1, 2, 3, 4};
Physical Line(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14};
Physical Line(2) = {12};
Physical Surface(0) = {1};