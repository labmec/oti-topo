/**
 * @file TPZHelmholtz.h
 * @brief Contains the TPZHelmholtz class which was firstly used as filter for topology optimization problems
 */

 #ifndef TPZHELMHOLTZ_H
 #define TPZHELMHOLTZ_H
 
 #include "TPZMatBase.h"
 #include "TPZMatSingleSpace.h"
 #include "TPZMatErrorSingleSpace.h"
 
 
 /**
  * @brief 
  */
 template<class TVar=STATE>
 class TPZHelmholtz :
     public TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>,TPZMatErrorSingleSpace<TVar>>{
     using TBase = TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>,TPZMatErrorSingleSpace<TVar>>;	
 public:
     //! Default constructor
     TPZHelmholtz() = default;
     /**
      * @brief Class constructor 
      * @param id material id
      * @param dim problem dimension
      * @param nstate number of state variables
      */
     TPZHelmholtz(int id, int dim, TVar r, int nstate=1);
 
     std::string Name() const override { return "TPZHelmholtz"; }
     
     /** @brief Solution indices of post-processing */
     enum ESolutionVars { ENone = 0, ESolution = 1 , EDerivative = 2};
     
     int Dimension() const  override { return this->fDim; }
     
     /** @brief Sets problem dimension */
     virtual void SetDimension(int dim) { this->fDim = dim; }
     
     int NStateVariables() const override { return this->fNStateVars; }
     void SetNStateVariables(int nstate) { this->fNStateVars = nstate; }
     void Contribute(const TPZMaterialDataT<TVar> &data,
                     REAL weight,
                     TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;
 
     void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                       TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                       TPZBndCondT<TVar> &bc) override;
     /** @brief To create another material of the same type */
     TPZMaterial * NewMaterial() const override;
     
     /** @brief It returns the variable index associated with the name */
     int VariableIndex(const std::string &name) const override;
     
     int NSolutionVariables(int var) const override;
 
     void Solution(const TPZMaterialDataT<TVar> &data,
                   int var, TPZVec<TVar> &solOut) override;
     
     void GetSolDimensions(uint64_t &u_len,
                           uint64_t &du_row,
                           uint64_t &du_col) const override;
 
     void Errors(const TPZMaterialDataT<TVar>&data,
                 TPZVec<REAL> &errors) override;
     
     virtual int ClassId() const override;
 protected:
     
     /** @brief Problem dimension */
     int fDim;
     
     /** @brief Number of state variables */
     int fNStateVars{1};
 
     /** @brief projection radius r */
     TVar fr2{0.0};
     
 };
 
 
 extern template class TPZHelmholtz<STATE>;
 extern template class TPZHelmholtz<CSTATE>;
 #endif
 