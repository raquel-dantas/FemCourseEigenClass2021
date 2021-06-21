//
//  TestOneDProblem.cpp
//  FemSC
//
//  Created by Eduardo Ferri on 08/17/15.
//
//
//TestOneDProblem cpp
/*
        Os testes foram preparados com um proposito educacional,
        recomenda-se que o aluno entenda a funcionalidade de cada
        teste e posteriormente use com seu código caso a caso
*/
//      Obs: O xmax e xmin estao tomados como 4 e 0, respectivamente,
//      caso estes valores sejam alterados, editar o teste TestNodes.
//
//
#include <iostream>
#include <math.h>
#include "GeoMesh.h"
#include "ReadGmsh.h"
#include "CompMesh.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "MathStatement.h"
#include "L2Projection.h"
#include "Analysis.h"
#include "IntRule.h"
#include "PostProcessTemplate.h"
#include "Poisson.h"

using std::cout;
using std::endl;
using std::cin;

void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv);

int main ()
{
    GeoMesh gmesh;
    ReadGmsh read;
    std::string filename("oneD.msh");
#ifdef MACOSX
    filename = "../"+filename;
#endif
    read.Read(gmesh,filename);

    CompMesh cmesh(&gmesh);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1,perm);
    mat1->SetDimension(1);
    
    auto force = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = 1.;
    };
    mat1->SetForceFunction(force);
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    L2Projection *bc_linha = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_point = new L2Projection(0,3,proj,val1,val2);
    std::vector<MathStatement *> mathvec = {0,mat1,bc_point,bc_linha};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(2);
    cmesh.AutoBuild();
    cmesh.Resequence();

    
    
    Analysis Analysis(&cmesh);
    Analysis.RunSimulation();
    
    PostProcessTemplate<Poisson> postprocess;
    postprocess.SetExact(exact);
    
    VecDouble errvec;
    errvec = Analysis.PostProcessError(std::cout, postprocess);
    
    
    return 0;
}
void exact(const VecDouble &point,VecDouble &val, MatrixDouble &deriv){

    deriv(0,0) = 4-point[0];
    val[0]=point[0]*(8.-point[0])/2.;
    return;
    double E=exp(1.0);
    VecDouble x(1);
    x[0]=point[0];
    
    val[0]=(30. + 100.*pow(E,100.) - 130.*pow(E,10.*x[0]) - 3*x[0] + 3*pow(E,100.)*x[0])/(10.*(-1. + pow(E,100.)));
    deriv(0,0)=(-3. + 3*pow(E,100) - 1300*pow(E,10*x[0]))/(10.*(-1 + pow(E,100)));
}


