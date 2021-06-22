/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Poisson.h"
///\cond
#include <functional>
#include <string.h>
///\endcond

Poisson::Poisson() {
}

Poisson::Poisson(int materialid, MatrixDouble &perm) {
    permeability = perm;
    this->SetMatID(materialid);
}

Poisson::Poisson(const Poisson &copy) {
    permeability = copy.permeability;
    forceFunction = copy.forceFunction;
}

Poisson &Poisson::operator=(const Poisson &copy) {
    permeability = copy.permeability;
    forceFunction = copy.forceFunction;
    return *this;
}

Poisson *Poisson::Clone() const {
    return new Poisson(*this);
}

Poisson::~Poisson() {
}

MatrixDouble Poisson::GetPermeability() const {
    return permeability;
}

void Poisson::SetPermeability(const MatrixDouble &perm) {
    permeability = perm;
}

int Poisson::NEvalErrors() const {
    return 3;
}

int Poisson::VariableIndex(const PostProcVar var) const {
    if (var == ENone) return ENone;
    if (var == ESol) return ESol;
    if (var == EDSol) return EDSol;
    if (var == EFlux) return EFlux;
    if (var == EForce) return EForce;
    if (var == ESolExact) return ESolExact;
    if (var == EDSolExact) return EDSolExact;
    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return -1;
}

Poisson::PostProcVar Poisson::VariableIndex(const std::string &name) {
    if (!strcmp("Sol", name.c_str())) return ESol;
    if (!strcmp("DSol", name.c_str())) return EDSol;
    if (!strcmp("Flux", name.c_str())) return EFlux;
    if (!strcmp("Force", name.c_str())) return EForce;
    if (!strcmp("SolExact", name.c_str())) return ESolExact;
    if (!strcmp("DSolExact", name.c_str())) return EDSolExact;
    else {
        std::cout << "variable not implemented" << std::endl;
    }
    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return ENone;
}

int Poisson::NSolutionVariables(const PostProcVar var) {
    if (var == ESol) return this->NState();
    if (var == EDSol) return this->Dimension();
    if (var == EFlux) return this->Dimension();
    if (var == EForce) return this->NState();
    if (var == ESolExact) return this->NState();
    if (var == EDSolExact) return this->Dimension();
    else {
        std::cout << "variable not implemented" << std::endl;
    }
    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return -1;
}

void Poisson::ContributeError(IntPointData &data, VecDouble &u_exact, MatrixDouble &du_exact, VecDouble &errors) const {
    errors.resize(NEvalErrors());
    errors.setZero();
    MatrixDouble gradu;
    MatrixDouble axes = data.axes;

    VecDouble u = data.solution;
    MatrixDouble dudx = data.dsoldx;

    this->Axes2XYZ(dudx, gradu, axes);

    double diff = 0.0;
    for (int i = 0; i < this->NState(); i++) {
        diff = (u[i] - u_exact[i]);
        errors[0] += diff*diff;
    }

    errors[1] = 0.;
    int dim = this->Dimension();
    int nstate = this->NState();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < nstate; j++) {
            diff = (gradu(i, j) - du_exact(i, j));
            errors[1] += diff*diff;
        }

    }
    errors[2] = errors[0] + errors[1];
}

void Poisson::Contribute(IntPointData &data, double weight, MatrixDouble &EK, MatrixDouble &EF) const {
    // Get data (shape functions and gradiaent of these shape functions) from the integration current point
    VecDouble& phi = data.phi;
    MatrixDouble& dphi = data.dphidx;
    // Permeability tensor from variational statement 
    MatrixDouble perm = this->GetPermeability();
    // Force function / Source function
    double f = 0.;
    auto force = this->GetForceFunction();
    if(force)
    {
        VecDouble resloc(1);
        force(data.x, resloc);
        f = resloc[0];
    }

    // Contribution on the element stiffiness matrix
    int nrows = EK.rows();
    int ncols = EK.cols();
    for(int i = 0; i < nrows; i++){
        // Gradient of shape functions
        MatrixDouble gradphi_i = dphi.col(i);
        // Flux = permeability * grad u
        MatrixDouble flux = perm*gradphi_i;
        for(int j = 0; j < ncols; j++){
            MatrixDouble gradphi_j = dphi.col(j);
            EK(i,j) += data.weight * data.detjac * Inner(flux, gradphi_j);
        }
    }
    // Contribution on the element load vector
    for(int i = 0; i < nrows; i++){
        EF(i,0) += data.weight * data.detjac * phi[i] * f;  
    } 

}

void Poisson::PostProcessSolution(const IntPointData &data, const int var, VecDouble &Solout) const {
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
            //+++++++++++++++++
            Solout.resize(nstate);
            for (int i = 0; i < nstate; i++) {
                Solout[i] = data.solution[i];
            }
            //+++++++++++++++++
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
        case 3: //EFlux
        {
            //+++++++++++++++++
            // Please implement me
            std::cout << "\nPLEASE IMPLEMENT ME\n" << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
            //+++++++++++++++++
        }
            break;

        case 4: //EForce
        {
            Solout.resize(nstate);
            VecDouble result(nstate);
            this->forceFunction(data.x, result);
            for (int i = 0; i < nstate; i++) {
                Solout[i] = result[i];
            }
        }
            break;

        case 5: //ESolExact
        {
            //+++++++++++++++++
            // Please implement me
            std::cout << "\nPLEASE IMPLEMENT ME\n" << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
            //+++++++++++++++++
        }
            break;
        case 6: //EDSolExact
        {
            //+++++++++++++++++
            // Please implement me
            std::cout << "\nPLEASE IMPLEMENT ME\n" << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
            //+++++++++++++++++
        }
            break;


        default:
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

double Poisson::Inner(MatrixDouble &S, MatrixDouble & T) const {
    double inner = 0;
    for (int i = 0; i < S.rows(); i++) {
        for (int j = 0; j < S.cols(); j++) {
            inner += S(i, j) * T(i, j);
        }
    }
    return inner;
}
