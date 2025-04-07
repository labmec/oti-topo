#include "TPZElasticity2DOtiTopoMP.h"
#include <pzaxestools.h>

TPZElasticity2DOtiTopoMP::TPZElasticity2DOtiTopoMP() : TPZRegisterClassId(&TPZElasticity2DOtiTopoMP::ClassId),
                                                       TPZMatCombinedSpacesT<STATE>(),
                                                       TPZElasticity2D() {}

TPZElasticity2DOtiTopoMP::TPZElasticity2DOtiTopoMP(int id, int dim) : TPZRegisterClassId(&TPZElasticity2DOtiTopoMP::ClassId),
                                                                      TPZElasticity2D(id) {}

TPZElasticity2DOtiTopoMP::TPZElasticity2DOtiTopoMP(int id, STATE E, STATE nu, STATE fx, STATE fy, int planestress) : TPZElasticity2D(id, E, nu, fx, fy, planestress) {
}

int TPZElasticity2DOtiTopoMP::VariableIndex(const std::string &name) const {
  if (!strcmp("topodensity", name.c_str())) return 100;
  return TPZElasticity2D::VariableIndex(name);
}

int TPZElasticity2DOtiTopoMP::NSolutionVariables(int var) const {
  switch (var) {
    case 100:
      return 1;
  }
  return TPZElasticity2D::NSolutionVariables(var);
}

int TPZElasticity2DOtiTopoMP::ClassId() const {
  return Hash("TPZElasticity2DOtiTopoMP") ^ TPZElasticity2D::ClassId() << 1;
}

TPZMaterial *TPZElasticity2DOtiTopoMP::NewMaterial() const {
  return new TPZElasticity2DOtiTopoMP(*this);
}

void TPZElasticity2DOtiTopoMP::Print(std::ostream &out) const {
  out << "Material Name: " << this->Name() << "\n";
  out << "Material Id: " << TPZElasticity2D::Id() << "\n";
  out << "Dimension: " << TPZElasticity2D::Dimension() << "\n\n";
}

void TPZElasticity2DOtiTopoMP::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                          REAL weight, TPZFMatrix<STATE> &ek,
                                          TPZFMatrix<STATE> &ef) {
  // Call the contribute method of father class TPZElasticity2D
  TPZFMatrix<STATE> temp_ek(ek.Rows(), ek.Cols());
  temp_ek.Zero();
  TPZElasticity2D::Contribute(datavec[0], weight, temp_ek, ef);

  // Get the value of the topo density material
  const STATE den = datavec[1].sol[0][0];

  // Multiply temp_ek by cte and add to ek
  const REAL degradation = den*den*den;
  temp_ek *= degradation;
  ek += temp_ek;
}

void TPZElasticity2DOtiTopoMP::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
  // Call the contributeBC of father class TPZElasticity2D
  TPZElasticity2D::ContributeBC(datavec[0], weight, ek, ef, bc);
}

void TPZElasticity2DOtiTopoMP::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                        int var, TPZVec<STATE> &sol) {
  if (var == 100) {
    sol[0] = datavec[1].sol[0][0];
    return;
  }
  TPZElasticity2D::Solution(datavec[0], var, sol);
}

void TPZElasticity2DOtiTopoMP::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
  DebugStop();  // not implemented
}

void TPZElasticity2DOtiTopoMP::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
  int nref = datavec.size();
  for (int i = 0; i < nref; i++) {
    datavec[i].SetAllRequirements(false);
    datavec[i].fNeedsNeighborSol = false;
    datavec[i].fNeedsNeighborCenter = false;
    datavec[i].fNeedsNormal = false;
    datavec[i].fNeedsHSize = false;
    datavec[i].fNeedsSol = true;
  }
}

void TPZElasticity2DOtiTopoMP::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
  // default is no specific data requirements
  int nref = datavec.size();
  for (int iref = 0; iref < nref; iref++) {
    datavec[iref].SetAllRequirements(false);
    datavec[iref].fNeedsSol = true;
  }
  datavec[0].fNeedsNormal = true;
  if (type == 50) {
    DebugStop();  // What is this?
    for (int iref = 0; iref < nref; iref++) {
      datavec[iref].fNeedsSol = false;
    }
  }
}