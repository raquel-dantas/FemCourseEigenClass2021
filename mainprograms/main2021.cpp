
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
#include "CompMesh.h"
#include "Poisson.h"
#include "NullStatement.h"
#include "L2Projection.h"
#include "CompElement.h"


using std::cout;
using std::endl;
using std::cin;

// Constant forcing function
void force(const VecDouble& x, VecDouble& result){
    result.resize(1);
    result[0] = 1;
}

int main (){
    // Setup of mesh reader
    ReadGmsh reader;
    // Geometric mesh
    GeoMesh mesh1;
    // Read mesh file from gmsh
    reader.Read(mesh1,"examples/mesh1.msh");
    // Print geometric mesh for verification
    VTKGeoMesh::PrintGMeshVTK(&mesh1,"mesh1.vtk");
    // Computational mesh created from geometrical mesh
    CompMesh cmesh(&mesh1);

    // Permeability tensor for Poisson problem
    MatrixDouble perm(2,2);
    perm.setIdentity();
    // Create a variational formulation
    Poisson *varform = new Poisson(1, perm);
    // Set variational formulation into mesh
    cmesh.SetMathStatement(1, varform);
    // Add forcing function to variational formulation
    varform->SetForceFunction(force);

    // Create computational elements and setup approximation space
    cmesh.AutoBuild();
    cmesh.Resequence();

    // Print computational mesh for verification
    VTKGeoMesh::PrintCMeshVTK(&cmesh, 2,"mesh1.vtk");


    // Compute element stiffness matrices
    int nelements = cmesh.GetElementVec().size();
    MatrixDouble ek;
    MatrixDouble ef;
    for(int i=0; i < nelements; i++){
        CompElement* ei = cmesh.GetElement(i);
        ei->CalcStiff(ek,ef);
        VecDouble t(4);
        t.setConstant(1);
        std::cout << '\n' << ek*t << '\n';    
    }

    
    return 0;
}
