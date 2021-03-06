H = 0.5;
h = 0.2;
//+
Point(1) = {0, 0, 0, H};
//+
Point(2) = {1, 0, 0, H};
//+
Point(3) = {1, 1, 0, H};
//+
Point(4) = {0, 1, 0, H};
//+
Point(5) = {0.5, 0.5, 0, h};
//+
Point(6) = {0.75, 0.5, 0, h};
//+
Point(7) = {0.5, 0.75, 0, h};
//+
Point(8) = {0.25, 0.5, 0, h};
//+
Point(9) = {0.5, 0.25, 0, h};
//+
Circle(1) = {6, 5, 7};
//+
Circle(2) = {7, 5, 8};
//+
Circle(3) = {8, 5, 9};
//+
Circle(4) = {9, 5, 6};
//+
Line(5) = {1, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 4};
//+
Line(8) = {4, 1};
//+
Line Loop(1) = {5, 6, 7, 8};
//+
Line Loop(2) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1, 2};
Physical Point(1) = {1, 2, 3, 4, 6, 7, 8, 9};
Physical Line(1) = {1, 2, 3, 4, 5, 7, 8};
Physical Line(2) = {6};
Physical Surface(0) = {1};
