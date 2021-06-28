/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Assemble.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "MathStatement.h"
#include "CompElement.h"



Assemble::Assemble() {
}

Assemble::Assemble(CompMesh *mesh) {
    cmesh = mesh;
}

Assemble::Assemble(const Assemble &copy) {
    cmesh = copy.cmesh;
}

Assemble &Assemble::operator=(const Assemble &copy) {
    cmesh = copy.cmesh;
    return *this;
}

void Assemble::SetMesh(CompMesh *mesh) {
    cmesh = mesh;
}

int64_t Assemble::NEquations() {
    int64_t neq = 0;
    int64_t i, ncon = cmesh->GetDOFVec().size();
    for (i = 0; i < ncon; i++) {
        DOF dof = cmesh->GetDOF(i);
        int64_t dofsize = dof.GetNShape() * dof.GetNState();
        neq += dofsize;
    }
    return neq;
}

void Assemble::OptimizeBandwidth() {    
}

void Assemble::Compute(MatrixDouble &globmat, MatrixDouble &rhs) {
    // Number of equations from the role mesh
    auto neq = NEquations();
    // Setup of the global matrix and rhs
    globmat.resize(neq, neq);
    globmat.setZero();
    rhs.resize(neq, 1);
    rhs.setZero();
    
    // Loop over computational elements
    int64_t nelem = cmesh->GetGeoMesh()->NumElements();
    for (int el = 0; el < nelem; el++) {
        CompElement *cel = cmesh->GetElement(el);
        // Number of shape functions defined over the element
        int nshape = cel->NShapeFunctions();
        // Number of state variables (1 for scalar, 3 for vectors)
        int nstate = cel->GetStatement()->NState();
        // Setup element stiffness matrix
        MatrixDouble ek(nstate * nshape, nstate * nshape);
        MatrixDouble ef(nstate * nshape, 1);
        ek.setZero();
        ef.setZero();
        //Compute element stiffness matrix and load vector
        cel->CalcStiff(ek, ef);

        // Loop over all DOFs
        int ndofs = cel->NDOF();
        for(int i = 0; i < ndofs; i++){
            // Get the DOF of i index
            DOF idof = cel->GetDOF(i);
            // Number of shape functions of the i index DOF (a DOF can have more than one shape functions)
            int nfun_i = idof.GetNShape();
            int64_t ifeq = idof.GetFirstEquation();
            // Loop over functions grouped in the DOF
            for(int ii = 0; ii < nfun_i*nstate; ii++){
                // Local index of the element matrix coefficient
                int iLocal = i + ii;
                // Global index associated to iLocal
                int iGlobal = ifeq + ii;
                
                // Add ef to the global rhs vector
                rhs(iGlobal,0) += ef(iLocal,0);

                // Loop over DOFs
                for(int j = 0; j < ndofs; j++){
                    // Get the DOF of j index
                    DOF jdof = cel->GetDOF(j);
                    // Number of shape functions of the j index DOF (a DOF can have more than one shape functions)
                    int nfun_j = jdof.GetNShape();
                    int64_t jfeq = jdof.GetFirstEquation();
                    // Loop over functions grouped in the DOF
                    for(int jj = 0; jj < nfun_j*nstate; jj++){
                        // Local index of the element matrix coefficient
                        int jLocal = j + jj;
                        // Global index associated to jLocal
                        int jGlobal = jfeq + jj;

                        // Add ek to the global stiffness matrix
                        globmat(iGlobal,jGlobal) += ek(iLocal,jLocal);

                    }
                }
            }
        }
        
    }
}

/**
        VecInt globalindex(ek.rows());

        int ndofs = cel->NDOF();
        for(int i = 0; i < ndofs; i++){
            DOF idof = cel->GetDOF(i);
            int nfun_i = idof.GetNShape();
            int64_t ifeq = idof.GetFirstEquation();

            for(int ii = 0; ii < nfun_i; ii++){
                int iLocal = i+ii;
                int iGlobal = ifeq + ii;
                globalindex[iLocal] = iGlobal;

            }
        }
        */