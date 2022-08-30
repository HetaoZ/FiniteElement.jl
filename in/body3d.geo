x0 = 0;
y0 = 0;
z0 = 0;
x1 = 1e-3;
y1 = 1e-3;
z1 = 1e-3;
lc = 5e-4;

Point(1) = {x0, y0, z0, lc};
Point(2) = {x1, y0, z0, lc};
Point(3) = {x0, y1, z0, lc};
Point(4) = {x1, y1, z0, lc};

Point(5) = {x0, y0, z1, lc};
Point(6) = {x1, y0, z1, lc};
Point(7) = {x0, y1, z1, lc};
Point(8) = {x1, y1, z1, lc};

Line(11) = {1,2};
Line(12) = {2,4};
Line(13) = {4,3};
Line(14) = {3,1};

Line(15) = {5,6};
Line(16) = {6,8};
Line(17) = {8,7};
Line(18) = {7,5};

Line(21) = {1,5};
Line(22) = {2,6};
Line(23) = {3,7};
Line(24) = {4,8};

Line Loop(9) = {5,6,7,8};
Plane Surface(10) = {9};
Recombine Surface(10);

Physical Line(101) = {7};
Physical Line(102) = {5};
Physical Line(103) = {8};
Physical Line(104) = {6};

Physical Surface("boundary") = {10};