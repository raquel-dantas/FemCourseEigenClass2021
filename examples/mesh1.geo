// Mesh IntroFem 1

Point(1) = {-1,-1,0};
Point(2) = {1,-1,0};
Point(3) = {1,1,0};
Point(4) = {-1,1,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};
Surface(1) = {1};

// To get structured triangular mesh 
Transfinite Curve{:} = 10;
Transfinite Surface{:};

//  To get structured quadilaterals
Recombine Surface{:};

Physical Surface("material",1) = Surface{:};