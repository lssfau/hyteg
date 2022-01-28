/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/solvers/Solver.hpp"

#include "PETScSparseMatrix.hpp"
#include "PETScVector.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

template < class OperatorType >
class PETScCGSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   PETScCGSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                  const uint_t&                              level,
                  const real_t                               relativeTolerance = 1e-30,
                  const real_t                               absoluteTolerance = 1e-12,
                  const PetscInt                             maxIterations     = std::numeric_limits< PetscInt >::max() )
   : allocatedLevel_( level )
   , petscCommunicator_( storage->getSplitCommunicatorByPrimitiveDistribution() )
   , num( "numerator", storage, level, level )
   , Amat( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ),
           numberOfGlobalDoFs< typename FunctionType::Tag >( *storage, level, petscCommunicator_ ),
           "Amat",
           petscCommunicator_ )
   , AmatNonEliminatedBC( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ),
                          numberOfGlobalDoFs< typename FunctionType::Tag >( *storage, level, petscCommunicator_ ),
                          "AmatNonEliminatedBC",
                          petscCommunicator_ )
   , xVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ), "xVec", petscCommunicator_ )
   , bVec( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ), "bVec", petscCommunicator_ )
   , nullspaceVec_( numberOfLocalDoFs< typename FunctionType::Tag >( *storage, level ), "nullspaceVec", petscCommunicator_ )
   , flag_( hyteg::All )
   , nullSpaceSet_( false )
   , reassembleMatrix_( false )
   {
      KSPCreate( petscCommunicator_, &ksp );
      KSPSetType( ksp, KSPCG );
      KSPSetTolerances( ksp, relativeTolerance, absoluteTolerance, PETSC_DEFAULT, maxIterations );
      KSPSetInitialGuessNonzero( ksp, PETSC_TRUE );
      KSPSetFromOptions( ksp );
   }

   ~PETScCGSolver()
   {
      KSPDestroy( &ksp );
      if ( nullSpaceSet_ )
         MatNullSpaceDestroy( &nullspace_ );
   }

   void reassembleMatrix( bool reassembleMatrix ) { reassembleMatrix_ = reassembleMatrix; }

   void setNullSpace( const FunctionType& nullspace )
   {
      nullSpaceSet_ = true;
      nullspaceVec_.createVectorFromFunction( nullspace, num, allocatedLevel_ );
      MatNullSpaceCreate( petscCommunicator_, PETSC_FALSE, 1, &nullspaceVec_.get(), &nullspace_ );
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level )
   {
      WALBERLA_CHECK_EQUAL( level, allocatedLevel_ );

      x.getStorage()->getTimingTree()->start( "PETSc CG Solver" );

      num.copyBoundaryConditionFromFunction( x );
      num.enumerate( level );

      xVec.createVectorFromFunction( x, num, level );
      bVec.createVectorFromFunction( b, num, level, All );

      if ( reassembleMatrix_ )
      {
         AmatNonEliminatedBC.zeroEntries();
         Amat.zeroEntries();
         AmatNonEliminatedBC.createMatrixFromOperator( A, level, num, All );
         Amat.createMatrixFromOperator( A, level, num, All );
      }
      else
      {
         AmatNonEliminatedBC.createMatrixFromOperatorOnce( A, level, num, All );
         Amat.createMatrixFromOperatorOnce( A, level, num, All );
      }
      MatCopy( AmatNonEliminatedBC.get(), Amat.get(), SAME_NONZERO_PATTERN );

      Amat.applyDirichletBCSymmetrically( x, num, bVec, level );
      if ( nullSpaceSet_ )
      {
         MatSetNullSpace( Amat.get(), nullspace_ );
      }
      KSPSetOperators( ksp, Amat.get(), Amat.get() );
      KSPGetPC( ksp, &pc );
      PCSetType( pc, PCHYPRE );

      KSPSolve( ksp, bVec.get(), xVec.get() );

      xVec.createFunctionFromVector( x, num, level, flag_ );

      x.getStorage()->getTimingTree()->stop( "PETSc CG Solver" );
   }

 private:
   uint_t                                                                                        allocatedLevel_;
   MPI_Comm                                                                                      petscCommunicator_;
   typename OperatorType::srcType::template FunctionType< idx_t >                                num;
   PETScSparseMatrix< OperatorType >                                                             Amat;
   PETScSparseMatrix< OperatorType >                                                             AmatNonEliminatedBC;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > xVec;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > bVec;
   PETScVector< typename FunctionType::valueType, OperatorType::srcType::template FunctionType > nullspaceVec_;

   KSP            ksp;
   PC             pc;
   MatNullSpace   nullspace_;
   hyteg::DoFType flag_;
   bool           nullSpaceSet_;
   bool           reassembleMatrix_;
};

} // namespace hyteg

#endif
