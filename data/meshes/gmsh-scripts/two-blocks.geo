// Setup a rectangular domain composed of two blocks

lc = 3;

// define six points, which will become the vertices of our two
// rectangular blocks
Point(1) = {0, 0, 0, lc};
Point(2) = {4, 0, 0, lc};
Point(3) = {9, 0, 0, lc};
Point(4) = {9, 3, 0, lc};
Point(5) = {4, 3, 0, lc};
Point(6) = {0, 3, 0, lc};

// construct edges of the blocks from the points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {5, 2};

// construct the two blocks from their four edges
Curve Loop(1) = {1, -7, 5, 6 };
Curve Loop(2) = {2, 3, 4, 7 };
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// assign names and physical tags to entities
Physical Point("Nodes on the Left", 1 ) = {1, 6};
Physical Point("Nodes on the Rigth", 3 ) = {3, 4};
Physical Point("Nodes in the Middle", 2 ) = {2, 5};

Physical Curve("Outer Boundary of Left Block", 1) = {1, 5, 6};
Physical Curve("Interface", 2) = {7};
Physical Curve("Outer Boundary of Right Block", 3) = {2, 3, 4};

Physical Surface( "Left Block", 4) = {1};
Physical Surface( "Right Block", 5) = {2};

// EOF
