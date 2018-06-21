#pragma once

#include <memory>

#include "PETScVector.hpp"
#include "PETScSparseMatrix.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"

#ifdef HHG_BUILD_WITH_PETSC

#if (PETSC_VERSION_MAJOR > 3) || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 9)
#define HHG_PCFactorSetMatSolverType PCFactorSetMatSolverType
#else
#define HHG_PCFactorSetMatSolverType PCFactorSetMatSolverPackage
#endif

namespace hhg{


template <typename ValueType, template <typename> class FunctionType, class OperatorType>
class PETScLUSolver {
public:

  PETScLUSolver(std::shared_ptr<FunctionType<PetscInt>> &numerator, uint_t localSize, uint_t globalSize)
      :num(numerator), Amat(localSize, globalSize), xVec(localSize), bVec(localSize)
  {
     for(uint_t i = num->getMinLevel(); i <= num->getMaxLevel(); ++i){
        hhg:communication::syncFunctionBetweenPrimitives(*num, i);
     }
     KSPCreate(PETSC_COMM_WORLD, &ksp);
  }

  ~PETScLUSolver(){
    KSPDestroy(&ksp);
  }

  void solve(OperatorType& A, FunctionType<ValueType>& x, FunctionType<ValueType>& b, FunctionType<ValueType>& r, size_t level, real_t tolerance, size_t maxiter, DoFType flag = All, bool printInfo = false) {

    bVec.createVectorFromFunction(b,*num.get(),level,All);

    if(Amat.createMatrixFromFunctionOnce(A, level, *num.get(), All))
    {
      Amat.applyDirichletBC(*num.get(),level);
      KSPSetOperators(ksp,Amat.get(),Amat.get());

      KSPGetPC(ksp,&pc);
      PCSetType(pc,PCLU);
      HHG_PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
      //PCFactorSetUpMatSolverPackage(pc); /* call MatGetFactor() to create F */
      //PCFactorGetMatrix(pc,&F);
    }


    //WALBERLA_LOG_INFO_ON_ROOT("Solving Linerar System")
    KSPSolve(ksp,bVec.get(),xVec.get());

    xVec.createFunctionFromVector(x,*num.get(),level,flag);

  }




private:
  std::shared_ptr<FunctionType<PetscInt>> num;
  PETScSparseMatrix<OperatorType,FunctionType> Amat;
  PETScVector<ValueType, FunctionType> xVec;
  PETScVector<ValueType, FunctionType> bVec;
  KSP ksp;
  PC pc;
  //Mat F; //factored Matrix



};



}

#endif
