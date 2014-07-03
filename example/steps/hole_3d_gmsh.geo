lc = 1;
r1 = 5;
r2 = 3;

Point(1) = {-r1,-r1,-r1,lc};
Point(2) = { r1,-r1,-r1,lc};
Point(3) = { r1, r1,-r1,lc};
Point(4) = {-r1, r1,-r1,lc};

Point(5) = {  0,  0,-r1,lc};

Point(6) = {  0,-r2,-r1,lc};
Point(7) = {  0, r2,-r1,lc};
Point(8) = { r2,  0,-r1,lc};
Point(9) = {-r2,  0,-r1,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {8, 5, 7};
Circle(6) = {7, 5, 9};
Circle(7) = {9, 5, 6};
Circle(8) = {6, 5, 8};

Line Loop(9)  = {1,2,3,4};
Line Loop(10) = {5,6,7,8};

Plane Surface(11) = {9, 10};
Recombine Surface {11};

Extrude {0,0,2*r1} {
  Surface{11}; Layers{2*r1/lc}; Recombine;
}
