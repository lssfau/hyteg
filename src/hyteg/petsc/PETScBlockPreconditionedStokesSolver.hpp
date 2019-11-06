#pragma once

#include <memory>

#include "hyteg/FunctionIterator.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "PETScSparseMatrix.hpp"
#include "PETScVector.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

template < class OperatorType >
class PETScBlockPreconditionedStokesSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType               FunctionType;
   typedef typename OperatorType::BlockPreconditioner_T BlockPreconditioner_T;

   PETScBlockPreconditionedStokesSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                                         const uint_t&                              level,
                                         const real_t                               tolerance = 1e-16,
                                         const PetscInt maxIterations = std::numeric_limits< PetscInt >::max() )
   : allocatedLevel_( level )
   , num( "numerator", storage, level, level )
   , Amat( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ),
           numberOfGlobalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , Pmat( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ),
           numberOfGlobalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , xVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , bVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , nullspaceVec_( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ) )
   , storage_( storage )
   , tolerance_( tolerance )
   , maxIterations_( maxIterations )
   , flag_( hyteg::All )
   , nullSpaceSet_( false )
   , blockPreconditioner_( storage, level, level )
   {}

   ~PETScBlockPreconditionedStokesSolver() {}

   void setNullSpace( const FunctionType& nullspace )
   {
      nullSpaceSet_ = true;
      nullspaceVec_.createVectorFromFunction( nullspace, num, allocatedLevel_ );
      MatNullSpaceCreate( walberla::MPIManager::instance()->comm(), PETSC_FALSE, 1, &nullspaceVec_.get(), &nullspace_ );
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level )
   {
      WALBERLA_CHECK_EQUAL( level, allocatedLevel_ );

      x.getStorage()->getTimingTree()->start( "PETSc block prec MinRes Solver" );

      x.getStorage()->getTimingTree()->start( "Setup" );

      x.getStorage()->getTimingTree()->start( "Index set setup" );
      num.enumerate( level );

      // gather index sets to split matrix into block matrix
      // therefore we need the row indices of the velocity and pressure
      std::vector< PetscInt > velocityIndices;
      std::vector< PetscInt > pressureIndices;

      gatherIndices( velocityIndices, pressureIndices, *storage_, level, num );

      std::sort( velocityIndices.begin(), velocityIndices.end() );
      std::sort( pressureIndices.begin(), pressureIndices.end() );
      ISCreateGeneral( walberla::mpi::MPIManager::instance()->comm(),
                       velocityIndices.size(),
                       velocityIndices.data(),
                       PETSC_COPY_VALUES,
                       &is_[0] );
      ISCreateGeneral( walberla::mpi::MPIManager::instance()->comm(),
                       pressureIndices.size(),
                       pressureIndices.data(),
                       PETSC_COPY_VALUES,
                       &is_[1] );

      x.getStorage()->getTimingTree()->stop( "Index set setup" );

      KSPCreate( walberla::MPIManager::instance()->comm(), &ksp );
      KSPSetType( ksp, KSPMINRES );
      KSPSetTolerances( ksp, tolerance_, tolerance_, PETSC_DEFAULT, maxIterations_ );
      // KSPSetInitialGuessNonzero( ksp, PETSC_TRUE );
      KSPSetFromOptions( ksp );

      x.getStorage()->getTimingTree()->start( "Vector copy" );
      xVec.createVectorFromFunction( x, num, level );
      bVec.createVectorFromFunction( b, num, level, All );
      x.getStorage()->getTimingTree()->stop( "Vector copy" );

      x.getStorage()->getTimingTree()->start( "Matrix assembly" );
      Amat.createMatrixFromFunctionOnce( A, level, num, All );
      Pmat.createMatrixFromFunctionOnce( blockPreconditioner_, level, num, All );
      x.getStorage()->getTimingTree()->stop( "Matrix assembly" );

      x.getStorage()->getTimingTree()->start( "Dirichlet BCs" );
      Amat.applyDirichletBCSymmetrically( x, num, bVec, level );
      Pmat.applyDirichletBCSymmetrically( num, level );
      x.getStorage()->getTimingTree()->stop( "Dirichlet BCs" );

      if ( nullSpaceSet_ )
      {
         MatSetNullSpace( Amat.get(), nullspace_ );
      }
      KSPSetOperators( ksp, Amat.get(), Pmat.get() );

      KSPGetPC( ksp, &pc );
      PCSetType( pc, PCFIELDSPLIT );
      PCFieldSplitSetType( pc, PC_COMPOSITE_ADDITIVE );
      PCFieldSplitSetIS( pc, "u", is_[0] );
      PCFieldSplitSetIS( pc, "p", is_[1] );

      PetscInt numSubKsps;
      PC       pc_u, pc_p;
      PCFieldSplitGetSubKSP( pc, &numSubKsps, &sub_ksps_ );

      KSPGetPC( sub_ksps_[0], &pc_u );
      KSPGetPC( sub_ksps_[1], &pc_p );

      // AMG
      PCSetType( pc_u, PCGAMG );
      PCGAMGSetType( pc_u, PCGAMGAGG );
      PCGAMGSetNSmooths( pc_u, 1 );

      // inv. lumped mass
      PCSetType( pc_p, PCJACOBI );

      x.getStorage()->getTimingTree()->stop( "Setup" );

      x.getStorage()->getTimingTree()->start( "Solve" );

      KSPSolve( ksp, bVec.get(), xVec.get() );

      x.getStorage()->getTimingTree()->stop( "Solve" );

      xVec.createFunctionFromVector( x, num, level, flag_ );

      x.getStorage()->getTimingTree()->stop( "PETSc block prec MinRes Solver" );

      PetscFree( sub_ksps_ );

      ISDestroy( &is_[0] );
      ISDestroy( &is_[1] );
      KSPDestroy( &ksp );
      if ( nullSpaceSet_ )
         MatNullSpaceDestroy( &nullspace_ );
   }

 private:
   void gatherIndices( std::vector< PetscInt >&            velocityIndices,
                       std::vector< PetscInt >&            pressureIndices,
                       const PrimitiveStorage&             storage,
                       uint_t                              level,
                       const P1StokesFunction< PetscInt >& numerator )
   {
      for ( auto dof : FunctionIterator< vertexdof::VertexDoFFunction< PetscInt > >( numerator.u, level ) )
      {
         velocityIndices.push_back( dof.value() );
      }

      for ( auto dof : FunctionIterator< vertexdof::VertexDoFFunction< PetscInt > >( numerator.v, level ) )
      {
         velocityIndices.push_back( dof.value() );
      }

      if ( storage.hasGlobalCells() )
      {
         for ( auto dof : FunctionIterator< vertexdof::VertexDoFFunction< PetscInt > >( numerator.w, level ) )
         {
            velocityIndices.push_back( dof.value() );
         }
      }
      for ( auto dof : FunctionIterator< vertexdof::VertexDoFFunction< PetscInt > >( numerator.p, level ) )
      {
         pressureIndices.push_back( dof.value() );
      }
   }

   void gatherIndices( std::vector< PetscInt >&                  velocityIndices,
                       std::vector< PetscInt >&                  pressureIndices,
                       const PrimitiveStorage&                   storage,
                       uint_t                                    level,
                       const P2P1TaylorHoodFunction< PetscInt >& numerator )
   {
      for ( auto dof : FunctionIterator< vertexdof::VertexDoFFunction< PetscInt > >( numerator.u.getVertexDoFFunction(), level ) )
      {
         velocityIndices.push_back( dof.value() );
      }
      for ( auto dof : FunctionIterator< EdgeDoFFunction< PetscInt > >( numerator.u.getEdgeDoFFunction(), level ) )
      {
         velocityIndices.push_back( dof.value() );
      }
      for ( auto dof : FunctionIterator< vertexdof::VertexDoFFunction< PetscInt > >( numerator.v.getVertexDoFFunction(), level ) )
      {
         velocityIndices.push_back( dof.value() );
      }
      for ( auto dof : FunctionIterator< EdgeDoFFunction< PetscInt > >( numerator.v.getEdgeDoFFunction(), level ) )
      {
         velocityIndices.push_back( dof.value() );
      }
      if ( storage.hasGlobalCells() )
      {
         for ( auto dof :
               FunctionIterator< vertexdof::VertexDoFFunction< PetscInt > >( numerator.w.getVertexDoFFunction(), level ) )
         {
            velocityIndices.push_back( dof.value() );
         }
         for ( auto dof : FunctionIterator< EdgeDoFFunction< PetscInt > >( numerator.w.getEdgeDoFFunction(), level ) )
         {
            velocityIndices.push_back( dof.value() );
         }
      }
      for ( auto dof : FunctionIterator< vertexdof::VertexDoFFunction< PetscInt > >( numerator.p, level ) )
      {
         pressureIndices.push_back( dof.value() );
      }
   }

   uint_t                                                                          allocatedLevel_;
   typename OperatorType::srcType::template FunctionType< PetscInt >               num;
   PETScSparseMatrix< OperatorType, OperatorType::srcType::template FunctionType > Amat;
   PETScSparseMatrix< BlockPreconditioner_T, BlockPreconditioner_T::srcType::template FunctionType > Pmat;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > xVec;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > bVec;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > nullspaceVec_;

   std::shared_ptr< PrimitiveStorage > storage_;

   real_t   tolerance_;
   PetscInt maxIterations_;

   BlockPreconditioner_T blockPreconditioner_;

   KSP          ksp;
   KSP*         sub_ksps_;
   PC           pc;
   IS           is_[2];
   MatNullSpace nullspace_;
   DoFType      flag_;
   bool         nullSpaceSet_;
};

} // namespace hyteg

#endif
