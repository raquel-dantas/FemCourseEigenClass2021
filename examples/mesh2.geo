// Extrusion mesh 

Point(1) = {-1,-1,0};
edge1[] = Extrude{2,0,0}{
    Point{1};
    Layers{10};
};

sur1[] = Extrude{0,2,0}{
    Line{1};
    Layers{10};
    Recombine;
};
