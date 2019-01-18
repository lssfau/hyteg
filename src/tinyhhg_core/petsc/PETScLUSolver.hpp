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
#if 0
  , inKernel( localSize )
#endif
  , flag_( hhg::All )
  {
     KSPCreate( walberla::MPIManager::instance()->comm(), &ksp );
  }

  ~PETScLUSolver(){
    KSPDestroy(&ksp);
  }

#if 0
  void setNullSpace( FunctionType & inKernel, const uint_t & level )
  {
    inKernel.createVectorFromFunction( inKernel, *num, level, All );
    VecNormalize(inKernel.get(), NULL);
    MatNullSpace nullspace;
    MatNullSpaceCreate( walberla::MPIManager::instance()->comm(), PETSC_FALSE, 1, &(inKernel.get()), &nullspace );
    MatSetNullSpace( Amat.get(), nullspace );
  }
#endif

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

    KSPSolve(ksp,bVec.get(),xVec.get());

    xVec.createFunctionFromVector(x,*num.get(),level,flag_);

  }




private:
  std::shared_ptr<typename OperatorType::srcType::template FunctionType<PetscInt>> num;
  PETScSparseMatrix<OperatorType,OperatorType::srcType::template FunctionType> Amat;
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> xVec;
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> bVec;
#if 0
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> inKernel;
#endif

  KSP ksp;
  PC pc;
  hhg::DoFType flag_;
  //Mat F; //factored Matrix



};



}

#endif
