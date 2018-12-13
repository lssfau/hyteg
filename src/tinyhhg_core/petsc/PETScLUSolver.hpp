#pragma once

#include <memory>

#include "PETScVector.hpp"
#include "PETScSparseMatrix.hpp"

#ifdef HHG_BUILD_WITH_PETSC

#if (PETSC_VERSION_MAJOR > 3) || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 9)
#define HHG_PCFactorSetMatSolverType PCFactorSetMatSolverType
#else
#define HHG_PCFactorSetMatSolverType PCFactorSetMatSolverPackage
#endif

namespace hhg{


template < class OperatorType >
class PETScLUSolver {
public:

  typedef typename OperatorType::srcType FunctionType;

  PETScLUSolver( std::shared_ptr< typename OperatorType::srcType::template FunctionType< PetscInt > >& numerator, uint_t localSize, uint_t globalSize )
  : num( numerator )
  , Amat( localSize, globalSize )
  , xVec( localSize )
  , bVec( localSize )
  , flag_( hhg::All )
  {
     KSPCreate( walberla::MPIManager::instance()->comm(), &ksp );
  }

  ~PETScLUSolver(){
    KSPDestroy(&ksp);
  }

  void solve(const OperatorType& A,const FunctionType& x,const FunctionType& b,const uint_t level) {

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

    xVec.createFunctionFromVector(x,*num.get(),level,flag_);

  }




private:
  std::shared_ptr<typename OperatorType::srcType::template FunctionType<PetscInt>> num;
  PETScSparseMatrix<OperatorType,OperatorType::srcType::template FunctionType> Amat;
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> xVec;
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> bVec;
  KSP ksp;
  PC pc;
  hhg::DoFType flag_;
  //Mat F; //factored Matrix



};



}

#endif
