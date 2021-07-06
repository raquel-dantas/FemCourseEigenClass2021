
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
    // double pi_y_over2 = M_PI*x[1]/2.;
    // result[0] = -2*std::cos(pi_y_over2)*std::sin(M_PI*x[1])+M_PI*M_PI/8.*(1+x[0])*(1+x[0])*(std::sin(pi_y_over2)+9*std::sin(3*pi_y_over2));
    // result[0] = 1./8.*(M_PI * M_PI)*(std::sin(M_PI * x[0]/2.)+9.*std::sin(3.*M_PI * x[0]/2.));
}
void f_sinsin(const VecDouble& x, VecDouble& result);

// Exact solution
void exactsol(const VecDouble& x, VecDouble& u, MatrixDouble& gradu);
void exactsol2(const VecDouble& x, VecDouble& u, MatrixDouble& gradu);
void sinsin(const VecDouble& x, VecDouble& u, MatrixDouble& gradu);

int main (){
    // Setup of mesh reader
    ReadGmsh reader;
    // Geometric mesh
    GeoMesh gmesh;
    // Read mesh file from gmsh
    reader.Read(gmesh,"examples/mesh_bc3.msh");
    // Print geometric mesh for verification
    VTKGeoMesh::PrintGMeshVTK(&gmesh,"gmesh_bc3.vtk");
    // Computational mesh created from geometrical mesh
    CompMesh cmesh(&gmesh);

    // Permeability tensor for Poisson problem
    MatrixDouble perm(2,2);
    perm.setIdentity();
    // Create a variational formulation
    Poisson *varform = new Poisson(0, perm);
    // Add forcing function to variational formulation
    varform->SetForceFunction(f_sinsin);
    // Create a BC
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
        proj.setZero();
        val1.setZero();
        val2.setZero();
        // val2(0,0) = 16;
    L2Projection *bc = new L2Projection(0,1,proj,val1,val2);
    bc->SetExactSolution(sinsin);
    // Set variational formulation into mesh
    std::vector<MathStatement*> ms {varform,bc};
    cmesh.SetMathVec(ms);

    // Create computational elements and setup approximation space
    cmesh.SetDefaultOrder(2);
    cmesh.AutoBuild();
    cmesh.Resequence();


    // Analysis 
    Analysis an(&cmesh);
    an.RunSimulation();

    // Print computational mesh for verification
    VTKGeoMesh::PrintCMeshVTK(&cmesh, 2,"result_sinsin.vtk");

    
    return 0;
}


void f_sinsin(const VecDouble& x, VecDouble& result){
    result.resize(1);
    result[0] = 2*(M_PI * M_PI) * std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
}

void exactsol(const VecDouble& x, VecDouble& u, MatrixDouble& gradu){
    u.resize(1);
    u[0] = std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
    // gradu.resize(2,1);
    // gradu(0,0) = M_PI * std::cos(M_PI * x[0]) * std::sin(M_PI * x[1]);
    // gradu(1,0) = M_PI * std::cos(M_PI * x[1]) * std::sin(M_PI * x[0]);
    u[0] = std::sin(M_PI * x[0]) * std::cos(M_PI * x[0]/2.);
}
void sinsin(const VecDouble& x, VecDouble& u, MatrixDouble& gradu){
    u.resize(1);
    u[0] = std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
    gradu.resize(2,1);
    gradu(0,0) = M_PI * std::cos(M_PI * x[0]) * std::sin(M_PI * x[1]);
    gradu(1,0) = M_PI * std::cos(M_PI * x[1]) * std::sin(M_PI * x[0]);
}

void exactsol2(const VecDouble& x, VecDouble& u, MatrixDouble& gradu){
    u.resize(1);
    u[0] = (x[0]+1)*(x[0]+1)*std::sin(M_PI*x[1])*std::cos(M_PI*x[1]/2.);
}