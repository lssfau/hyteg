#pragma once

#include <memory>

#include "tinyhhg_core/solvers/Solver.hpp"

#include "PETScSparseMatrix.hpp"
#include "PETScVector.hpp"

#ifdef HHG_BUILD_WITH_PETSC

namespace hhg {

template < class OperatorType >
class PETScMinResSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   PETScMinResSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                      const uint_t&                              level,
                      const real_t                               tolerance     = 1e-16,
                      const PetscInt                             maxIterations = std::numeric_limits< PetscInt >::max() )
   : allocatedLevel_( level )
   , num( "numerator", storage, level, level )
   , Amat( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ),
           numberOfGlobalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , AmatNonEliminatedBC( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ),
                          numberOfGlobalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , xVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , bVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , nullspaceVec_( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , flag_( hhg::All )
   , nullSpaceSet_( false )
   {
      num.enumerate( level );
      KSPCreate( walberla::MPIManager::instance()->comm(), &ksp );
      KSPSetType( ksp, KSPMINRES );
      KSPSetTolerances( ksp, tolerance, tolerance, PETSC_DEFAULT, maxIterations );
      KSPSetInitialGuessNonzero( ksp, PETSC_TRUE );
   }

   ~PETScMinResSolver()
   {
      KSPDestroy( &ksp );
      if ( nullSpaceSet_ )
         MatNullSpaceDestroy( &nullspace_ );
   }

   void setNullSpace( const FunctionType& nullspace )
   {
      nullSpaceSet_ = true;
      nullspaceVec_.createVectorFromFunction( nullspace, num, allocatedLevel_ );
      MatNullSpaceCreate( walberla::MPIManager::instance()->comm(), PETSC_FALSE, 1, &nullspaceVec_.get(), &nullspace_ );
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level )
   {
      WALBERLA_CHECK_EQUAL( level, allocatedLevel_ );

      x.getStorage()->getTimingTree()->start( "PETSc MinRes Solver" );

      xVec.createVectorFromFunction( x, num, level );
      bVec.createVectorFromFunction( b, num, level, All );

      AmatNonEliminatedBC.createMatrixFromFunctionOnce( A, level, num, All );
      Amat.createMatrixFromFunctionOnce( A, level, num, All );
      MatCopy( AmatNonEliminatedBC.get(), Amat.get(), DIFFERENT_NONZERO_PATTERN );

      Amat.applyDirichletBCSymmetrically( x, num, bVec, level );
      if ( nullSpaceSet_ )
      {
         MatSetNullSpace( Amat.get(), nullspace_ );
      }
      KSPSetOperators( ksp, Amat.get(), Amat.get() );
      KSPGetPC( ksp, &pc );
      PCSetType( pc, PCNONE );

      KSPSolve( ksp, bVec.get(), xVec.get() );

      xVec.createFunctionFromVector( x, num, level, flag_ );

      x.getStorage()->getTimingTree()->stop( "PETSc MinRes Solver" );
   }

 private:
   uint_t                                                                                        allocatedLevel_;
   typename OperatorType::srcType::template FunctionType< PetscInt >                             num;
   PETScSparseMatrix< OperatorType, OperatorType::srcType::template FunctionType >               Amat;
   PETScSparseMatrix< OperatorType, OperatorType::srcType::template FunctionType >               AmatNonEliminatedBC;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > xVec;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > bVec;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > nullspaceVec_;

   KSP          ksp;
   PC           pc;
   MatNullSpace nullspace_;
   hhg::DoFType flag_;
   //Mat F; //factored Matrix
   bool nullSpaceSet_;
};

} // namespace hhg

#endif
