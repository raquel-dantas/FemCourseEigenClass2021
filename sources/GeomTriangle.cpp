/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTriangle.h"

GeomTriangle::GeomTriangle() {
}

GeomTriangle::~GeomTriangle() {
}

GeomTriangle::GeomTriangle(const GeomTriangle &copy) {
    fNodeIndices = copy.fNodeIndices;

}

GeomTriangle& GeomTriangle::operator=(const GeomTriangle& copy) {
    fNodeIndices = copy.fNodeIndices;

    return *this;
}

void GeomTriangle::Shape(const VecDouble& xi, VecDouble& phi, MatrixDouble& dphi) {
    if(xi.size() != Dimension ) DebugStop();
    phi.resize(nCorners);
    phi[0] = 1. - xi[0] - xi[1];
    phi[1] = xi[0];
    phi[2] = xi[1];

     // dphi(i,j) represents the ith derivative of function j
    dphi.resize(Dimension,nCorners);
    dphi(0,0) = -1.;
    dphi(1,0) = -1.;
    dphi(0,1) = 1.;
    dphi(1,1) = 0.;
    dphi(0,2) = 0.;
    dphi(1,2) = 1.;
}

void GeomTriangle::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    
    VecDouble phi;
    MatrixDouble dphi;
    Shape(xi, phi, dphi);
    int nnodes = NumNodes();
    int dim = NodeCo.rows();
    // x.resize(dim);
    x.setZero();

    // NodeCo is a matrix dim x Nnodes (3 by 2 for a line). Taking NodeCo(i,j) gives the ith coordinate of the jth node

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < nnodes; j++)
        {
            x[i] += phi[j] * NodeCo(i, j);
        }
    }
}

void GeomTriangle::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();

    VecDouble phi;
    MatrixDouble dphi;
    Shape(xi, phi, dphi);
    int nnodes = NumNodes();
    int masterdim = Dimension;
    int compdim = NodeCo.rows();
    // x.resize(compdim);
    x.setZero();
    gradx.resize(compdim,masterdim);
    gradx.setZero();

    // NodeCo is a matrix dim x Nnodes (3 by 2 for a line). Taking NodeCo(i,j) gives the ith coordinate of the jth node
    for (int k = 0; k < nnodes; k++)
    {
        for (int i = 0; i < compdim; i++)
        {
            x[i] += phi[k] * NodeCo(i, k);
            for (int j = 0; j < masterdim; j++)
            {
                gradx(i, j) += NodeCo(i, k) * dphi(j, k);
            }
        }
    }
}

void GeomTriangle::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) DebugStop();
    fNodeIndices = nodes;
}

void GeomTriangle::GetNodes(VecInt &nodes) const  {
    nodes = fNodeIndices;
}

int GeomTriangle::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomTriangle::NumNodes() {
    return nCorners;
}

GeoElementSide GeomTriangle::Neighbour(int side)  const {
    return fNeighbours[side];
}

void GeomTriangle::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
