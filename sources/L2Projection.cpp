/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "L2Projection.h"
///\cond
#include <string.h>
///\endcond

L2Projection::L2Projection() {
}

L2Projection::L2Projection(int bctype, int materialid, MatrixDouble &proj, MatrixDouble Val1, MatrixDouble Val2) {
    projection = proj;
    BCType = bctype;
    BCVal1 = Val1;
    BCVal2 = Val2;
    this->SetMatID(materialid);
}

L2Projection::L2Projection(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;
    BCVal1 = copy.BCVal1;
    BCVal2 = copy.BCVal2;

}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;
    BCVal1 = copy.BCVal1;
    BCVal2 = copy.BCVal2;
    return *this;
}

L2Projection *L2Projection::Clone() const {
    return new L2Projection(*this);
}

L2Projection::~L2Projection() {
}

MatrixDouble L2Projection::GetProjectionMatrix() const {
    return projection;
}

void L2Projection::SetProjectionMatrix(const MatrixDouble &proj) {
    projection = proj;
}

void L2Projection::Contribute(IntPointData &data, double weight, MatrixDouble &EK, MatrixDouble &EF) const {
    // Get number of state variables
    int nstate = this->NState();
    if(nstate != 1)
    {
        std::cout << "Please implement me\n";
        DebugStop();
    }
    // Get number of shape functions
    auto nshape = data.phi.size();
    if(EK.rows() != nstate*nshape || EF.rows() != nstate*nshape)
    {
        DebugStop();
    }
    // Get the value of the user defined BC
    VecDouble bc_value(nstate);
    bc_value[0] = Val2()(0,0);

    // If exact solution is known, BCs are whatever we compute from it
    MatrixDouble gradu_exact(data.x.size(), nstate);
    gradu_exact.setZero();
    if(SolutionExact)
    {
        SolutionExact(data.x, bc_value, gradu_exact);
        // if(GetBCType() == EBCType::Neumann){
        //     MatrixDouble flux_exact = perm * gradu_exact;
        //     double normal_flux = Inner(flux_exact,normal_vector);
        //     bc_value[0] = normal_flux;
        // }
    }
    // Get data (shape functions) from the integration point
    VecDouble& phi = data.phi;

    switch (this->GetBCType()) {

        case 0://Dirichlet
        {
            int nrows = EK.rows();
            int ncols = EK.cols();
            // BCs are taken in the form of boundary integrals
            for(int i = 0; i < nrows; i++){
                // Bignumber because of panalty method
                EF(i,0) += gBigNumber * data.weight * data.detjac * phi[i] * bc_value[0];
                for(int j = 0; j < ncols; j++){
                    EK(i,j) += gBigNumber * data.weight * data.detjac * phi[i] * phi[j];
                }
            }
            break;
        }

        case 1://Neumann
        {
            // BCs are taken in the form of boundary integrals
            int nrows = EK.rows();
            for(int i = 0; i < nrows; i++){
                EF(i,0) += data.weight * data.detjac * phi[i] * bc_value[0];  
            } 
            break;
        }

        default:
        {
            std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
        }
    }
    //+++++++++++++++++
}

int L2Projection::NEvalErrors() const {
    return 3;
}

void L2Projection::ContributeError(IntPointData &data, VecDouble &u_exact, MatrixDouble &du_exact, VecDouble &errors) const {
    return;
}

int L2Projection::VariableIndex(const PostProcVar var) const {
    if (var == ESol) return ESol;
    if (var == EDSol) return EDSol;

    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return 0;
}

L2Projection::PostProcVar L2Projection::VariableIndex(const std::string & name) {
    if (!strcmp("Solution", name.c_str())) return ESol;
    if (!strcmp("Derivative", name.c_str())) return EDSol;

    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return ENone;
}

int L2Projection::NSolutionVariables(const PostProcVar var) {
    if (var == ESol) return this->NState();
    if (var == EDSol) return this->NState();

    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return 0;
}

void L2Projection::PostProcessSolution(const IntPointData &data, const int var, VecDouble &Solout) const {
    VecDouble sol = data.solution;
    int solsize = sol.size();
    int rows = data.dsoldx.rows();
    int cols = data.dsoldx.cols();
    MatrixDouble gradu(rows, cols);
    gradu = data.dsoldx;

    int nstate = this->NState();

    switch (var) {
        case 0: //None
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();
        }

        case 1: //ESol
        {
            Solout.resize(nstate);
            for (int i = 0; i < nstate; i++) {
                Solout[i] = data.solution[i];
            }
        }
            break;

        case 2: //EDSol
        {
           Solout.resize(rows*cols);
            for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                Solout[cols*i+j] = data.dsoldx.coeff(i,j);
            }
            }
        }
            break;

        default:
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}
