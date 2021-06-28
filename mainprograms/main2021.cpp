
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
    // result[0] = (M_PI * M_PI)/2. * std::cos(M_PI * x[0]/2.) * std::cos(M_PI * x[1]/2.);
    result[0] = 2*(M_PI * M_PI) * std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
}
// Exact solution
void exactsol(const VecDouble& x, VecDouble& u, MatrixDouble& gradu);

int main (){
    // Setup of mesh reader
    ReadGmsh reader;
    // Geometric mesh
    GeoMesh gmesh;
    // Read mesh file from gmsh
    reader.Read(gmesh,"examples/mesh_bc3.msh");
    // Print geometric mesh for verification
    VTKGeoMesh::PrintGMeshVTK(&gmesh,"gmesh_bc.vtk");
    // Computational mesh created from geometrical mesh
    CompMesh cmesh(&gmesh);

    // Permeability tensor for Poisson problem
    MatrixDouble perm(2,2);
    perm.setIdentity();
    // Create a variational formulation
    Poisson *varform = new Poisson(0, perm);
    // Add forcing function to variational formulation
    varform->SetForceFunction(force);
    // Create a BC
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
        proj.setZero();
        val1.setZero();
        val2.setZero();
        val2(0,0) = 16;
    L2Projection *bc = new L2Projection(0,1,proj,val1,val2);
    bc->SetExactSolution(exactsol);
    // Set variational formulation into mesh
    std::vector<MathStatement*> ms {varform,bc};
    cmesh.SetMathVec(ms);

    // Create computational elements and setup approximation space
    cmesh.AutoBuild();
    cmesh.Resequence();


    // Analysis 
    Analysis an(&cmesh);
    an.RunSimulation();

    // Print computational mesh for verification
    VTKGeoMesh::PrintCMeshVTK(&cmesh, 2,"result.vtk");

    return 0;
}




void exactsol(const VecDouble& x, VecDouble& u, MatrixDouble& gradu){
    u.resize(1);
    u[0] = std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
    gradu.resize(2,1);
    gradu(0,0) = M_PI * std::cos(M_PI * x[0]) * std::sin(M_PI * x[1]);
    gradu(0,0) = M_PI * std::cos(M_PI * x[1]) * std::sin(M_PI * x[0]);
}