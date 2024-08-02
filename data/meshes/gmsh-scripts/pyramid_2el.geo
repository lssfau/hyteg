// Generate a simple mesh of a straight pyramid consisting of two tetrahedrons

lc = 10;

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {1.0, 0.0, 0.0, lc};
Point(3) = {0.0, 1.0, 0.0, lc};
Point(4) = {0.5, 0.5, 1.0, lc};
Point(5) = {1.0, 1.0, 0.0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};
Line(4) = {1,4};
Line(5) = {2,4};
Line(6) = {3,4};
Line(7) = {2,5};
Line(8) = {5,3};
Line(9) = {5,4};

Curve Loop(1) = {1,2,3};
Curve Loop(2) = {1,5,-4};
Curve Loop(3) = {2,6,-5};
Curve Loop(4) = {-3,6,-4};
Curve Loop(5) = {7,8,-2};
Curve Loop(6) = {7,9,-5};
Curve Loop(7) = {9,-6,-8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};

Surface Loop(1) = {1,2,3,4};
Surface Loop(2) = {3,5,6,7};
Volume(1) = {1};
Volume(2) = {2};

Physical Point(1) = {1,2,3,4,5};
Physical Curve(1) = {1,2,3,4,5,6,7,8,9};
Physical Surface(1) = {1,2,3,4,5,6,7};
Physical Volume(0) = {1,2};

// EOF


