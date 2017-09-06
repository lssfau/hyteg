#pragma once

#include <memory>

#include "PETScVector.hpp"
#include "PETScSparseMatrix.hpp"

namespace hhg{


template <class Functiontype,class Operatortype>
class PETScLUSolver {
public:

  PETScLUSolver(std::shared_ptr<Functiontype> &numerator,uint_t size)
      :num(numerator),Amat(size),vec(size)
  {
    KSPCreate(PETSC_COMM_WORLD, &ksp);
  }

  ~PETScLUSolver(){
    KSPDestroy(&ksp);
  }

  void solve(Operatortype& A, Functiontype& x, Functiontype& b, Functiontype& r, size_t level, real_t tolerance, size_t maxiter, DoFType flag = All, bool printInfo = false) {

    vec.createVectorFromFunction(b,*num.get(),level,flag);

    if(Amat.createMatrixFromFunctionOnce(A, level, *num.get(), flag))
    {
      Amat.applyDirichletBC(*num.get(),level);
      KSPSetOperators(ksp,Amat.get(),Amat.get());

      KSPGetPC(ksp,&pc);
      PCSetType(pc,PCLU);
      PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
      //PCFactorSetUpMatSolverPackage(pc); /* call MatGetFactor() to create F */
      //PCFactorGetMatrix(pc,&F);
    }


    //WALBERLA_LOG_INFO_ON_ROOT("Solving Linerar System")
    KSPSolve(ksp,vec.get(),vec.get());

    vec.createFunctionFromVector(x,*num.get(),level,flag);

  }




private:
  std::shared_ptr<Functiontype> num;
  PETScSparseMatrix<Operatortype,Functiontype> Amat;
  PETScVector<Functiontype> vec;
  KSP ksp;
  PC pc;
  //Mat F; //factored Matrix



};



}
