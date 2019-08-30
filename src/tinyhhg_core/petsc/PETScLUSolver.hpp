#pragma once

#include <memory>

#include "PETScVector.hpp"
#include "PETScSparseMatrix.hpp"
#include "tinyhhg_core/solvers/Solver.hpp"

#ifdef HHG_BUILD_WITH_PETSC

#if (PETSC_VERSION_MAJOR > 3) || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 9)
#define HHG_PCFactorSetMatSolverType PCFactorSetMatSolverType
#else
#define HHG_PCFactorSetMatSolverType PCFactorSetMatSolverPackage
#endif

namespace hyteg {


template < class OperatorType >
class PETScLUSolver : public Solver< OperatorType >
{
public:

  typedef typename OperatorType::srcType FunctionType;

  PETScLUSolver( const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & level )
  : allocatedLevel_( level )
  , num( "numerator", storage, level, level )
  , Amat( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ), numberOfGlobalDoFs< typename FunctionType::Tag >( *storage, level ) )
  , xVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
  , bVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
#if 0
  , inKernel( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
#endif
  , flag_( hyteg::All )
  {
    num.enumerate( level );
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

  void setConstantNullSpace()
  {
    MatNullSpace nullspace;
    MatNullSpaceCreate( walberla::MPIManager::instance()->comm(), PETSC_TRUE, 0, NULL, &nullspace );
    MatSetNullSpace( Amat.get(), nullspace );
  }

  void solve(const OperatorType& A,const FunctionType& x, const FunctionType& b, const uint_t level)
  {
    WALBERLA_CHECK_EQUAL( level, allocatedLevel_ );

    x.getStorage()->getTimingTree()->start( "PETSc LU Solver" );

    b.assign({1.0}, {x}, level, DirichletBoundary);

    bVec.createVectorFromFunction(b,num,level,All);

    x.getStorage()->getTimingTree()->start( "Matrix assembly and solver setup" );
    const bool matrixAssembledForTheFirstTime = Amat.createMatrixFromFunctionOnce(A, level, num, All);

    if( matrixAssembledForTheFirstTime )
    {
      Amat.applyDirichletBC(num,level);
      KSPSetOperators(ksp,Amat.get(),Amat.get());

      KSPGetPC(ksp,&pc);
      PCSetType(pc,PCLU);
      HHG_PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
      //PCFactorSetUpMatSolverPackage(pc); /* call MatGetFactor() to create F */
      //PCFactorGetMatrix(pc,&F);
    }
    x.getStorage()->getTimingTree()->stop( "Matrix assembly and solver setup" );

    x.getStorage()->getTimingTree()->start( "Solver" );
    KSPSolve(ksp,bVec.get(),xVec.get());
    x.getStorage()->getTimingTree()->stop( "Solver" );

    xVec.createFunctionFromVector(x,num,level,flag_);

    x.getStorage()->getTimingTree()->stop( "PETSc LU Solver" );
  }




private:
  uint_t allocatedLevel_;
  typename OperatorType::srcType::template FunctionType<PetscInt> num;
  PETScSparseMatrix<OperatorType,OperatorType::srcType::template FunctionType> Amat;
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> xVec;
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> bVec;
#if 0
  PETScVector<typename FunctionType::valueType, OperatorType::srcType::template FunctionType> inKernel;
#endif

  KSP ksp;
  PC pc;
  hyteg::DoFType flag_;
  //Mat F; //factored Matrix



};



}

#endif
