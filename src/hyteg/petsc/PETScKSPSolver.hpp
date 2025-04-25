/*
 * Copyright (c) 2017-2025 Andreas Burkhart, Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <memory>

#include "hyteg/petsc/PETScSolverOptions.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "PETScSparseMatrix.hpp"
#include "PETScVector.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

PetscErrorCode PETScKSPSolverRemoveNullSpace( MatNullSpace nullspace, Vec vec, void* fct )
{
   typedef std::function< void( Vec ) > fType;
   typedef fType*                       fTypePtr;
   ( *( reinterpret_cast< fTypePtr >( fct ) ) )( vec );
   return PETSC_SUCCESS;
}

template < class OperatorType >
class PETScKSPSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType                                 FunctionType;
   typedef typename OperatorType::srcType::template FunctionType< idx_t > FunctionIdxType;

   PETScKSPSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                   const uint_t&                              level,
                   const PETScSolverOptions&                  solverOptions,
                   const real_t                               relativeTolerance      = 1e-30,
                   const real_t                               absoluteTolerance      = 1e-12,
                   const PetscInt                             maxIterations          = std::numeric_limits< PetscInt >::max(),
                   bool                                       alwaysReassembleMatrix = false,
                   bool                                       zeroInitialGuess       = false,
                   const std::string                          prefix                 = "" )
   : PETScKSPSolver( storage,
                     level,
                     solverOptions,
                     FunctionIdxType( "numerator", storage, level, level ),
                     relativeTolerance,
                     absoluteTolerance,
                     maxIterations,
                     alwaysReassembleMatrix,
                     zeroInitialGuess,
                     prefix )
   {}

   PETScKSPSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                   const uint_t&                              level,
                   const PETScSolverOptions&                  solverOptions,
                   const FunctionIdxType&                     numerator,
                   const real_t                               relativeTolerance      = 1e-30,
                   const real_t                               absoluteTolerance      = 1e-12,
                   const PetscInt                             maxIterations          = std::numeric_limits< PetscInt >::max(),
                   bool                                       alwaysReassembleMatrix = false,
                   bool                                       zeroInitialGuess       = false,
                   const std::string                          prefix                 = "" )
   : allocatedLevel_( level )
   , petscCommunicator_( storage->getSplitCommunicatorByPrimitiveDistribution() )
   , num( numerator )
   , Amat( "Amat", petscCommunicator_ )
   , AmatNonEliminatedBC( "AmatNonEliminatedBC", petscCommunicator_ )
   , xVec( "xVec", petscCommunicator_ )
   , bVec( "bVec", petscCommunicator_ )
   , nullspaceVec_( "nullspaceVec", petscCommunicator_ )
   , flag_( hyteg::All )
   , nullSpaceSet_( false )
   , alwaysReassembleMatrix_( alwaysReassembleMatrix )
   , disableApplicationBC_( false )
   , firstAssemble_( true )
   {
      num.enumerate( level );
      KSPCreate( petscCommunicator_, &ksp );
      if ( prefix != "" )
      {
         KSPSetOptionsPrefix( ksp, prefix.c_str() );
      }

      KSPSetTolerances( ksp, relativeTolerance, absoluteTolerance, PETSC_DEFAULT, maxIterations );
      if ( zeroInitialGuess )
      {
         KSPSetInitialGuessNonzero( ksp, PETSC_FALSE );
      }
      else
      {
         KSPSetInitialGuessNonzero( ksp, PETSC_TRUE );
      }

      KSPSetType( ksp, solverOptions.getKspType().c_str() );
      KSPGetPC( ksp, &pc );
      PCSetType( pc, solverOptions.getPcType().c_str() );
      solverOptions.applyOptions( ksp, prefix );
      KSPSetFromOptions( ksp );
   }

   ~PETScKSPSolver()
   {
      KSPDestroy( &ksp );
      if ( nullSpaceSet_ )
         MatNullSpaceDestroy( &nullspace_ );
   }

   void setAlwaysReassembleMatrix( bool alwaysReassembleMatrix ) { alwaysReassembleMatrix_ = alwaysReassembleMatrix; }

   void disableApplicationBC( bool disableApplicationBC ) { disableApplicationBC_ = disableApplicationBC; }

   void reassembleMatrix( const OperatorType& A, const uint_t level )
   {
      AmatNonEliminatedBC.zeroEntries();
      Amat.zeroEntries();
      AmatNonEliminatedBC.createMatrixFromOperator( A, level, num, All );
      Amat.createMatrixFromOperator( A, level, num, All );

      firstAssemble_ = false;
   }

   void setNullSpace( const FunctionType& nullspace )
   {
      nullSpaceSet_ = true;
      nullspaceVec_.createVectorFromFunction( nullspace, num, allocatedLevel_ );
      real_t norm = 0;
      VecNormalize( nullspaceVec_.get(), &norm );
      MatNullSpaceCreate( petscCommunicator_, PETSC_FALSE, 1, &nullspaceVec_.get(), &nullspace_ );
   }

   MatNullSpace*    getNullSpace() { return &nullspace_; }
   FunctionIdxType& getEnumerator() { return num; }

   void setNullSpaceFunction( std::shared_ptr< std::function< void( Vec ) > > nullSpaceFunction )
   {
      nullSpaceFunction_ = nullSpaceFunction;
      MatNullSpaceCreate( petscCommunicator_, PETSC_FALSE, 0, nullptr, &nullspace_ );
      MatNullSpaceSetFunction( nullspace_, hyteg::PETScKSPSolverRemoveNullSpace, nullSpaceFunction_.get() );
      nullSpaceSet_ = true;
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level )
   {
      WALBERLA_CHECK_EQUAL( level, allocatedLevel_ );

      x.getStorage()->getTimingTree()->start( "PETSc KSP Solver" );

      num.copyBoundaryConditionFromFunction( x );

      xVec.createVectorFromFunction( x, num, level );
      bVec.createVectorFromFunction( b, num, level, All );

      if ( firstAssemble_ || alwaysReassembleMatrix_ )
      {
         reassembleMatrix( A, level );
      }

      MatCopy( AmatNonEliminatedBC.get(), Amat.get(), DIFFERENT_NONZERO_PATTERN );

      if ( !disableApplicationBC_ )
         Amat.applyDirichletBCSymmetrically( x, num, bVec, level );
      if ( nullSpaceSet_ )
      {
         MatSetNullSpace( Amat.get(), nullspace_ );
      }
      KSPSetOperators( ksp, Amat.get(), Amat.get() );

      KSPSolve( ksp, bVec.get(), xVec.get() );

      xVec.createFunctionFromVector( x, num, level, flag_ );

      x.getStorage()->getTimingTree()->stop( "PETSc KSP Solver" );
   }

 private:
   uint_t                                                                                        allocatedLevel_;
   MPI_Comm                                                                                      petscCommunicator_;
   FunctionIdxType                                                                               num;
   PETScSparseMatrix< OperatorType >                                                             Amat;
   PETScSparseMatrix< OperatorType >                                                             AmatNonEliminatedBC;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > xVec;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > bVec;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > nullspaceVec_;
   std::shared_ptr< std::function< void( Vec ) > >                                               nullSpaceFunction_;

   KSP            ksp;
   PC             pc;
   MatNullSpace   nullspace_;
   hyteg::DoFType flag_;
   bool           nullSpaceSet_;
   bool           alwaysReassembleMatrix_;
   bool           disableApplicationBC_;
   bool           firstAssemble_;
};

} // namespace hyteg

#endif
