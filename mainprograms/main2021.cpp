
///\cond
#include <iostream>
#include <math.h>
///\endcond
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"

using std::cout;
using std::endl;
using std::cin;

int main (){
    ReadGmsh reader;
    GeoMesh mesh1;
    reader.Read(mesh1,"examples/mesh1.msh");
    VTKGeoMesh::PrintGMeshVTK(&mesh1,"mesh1.vtk");    
    return 0;
}
