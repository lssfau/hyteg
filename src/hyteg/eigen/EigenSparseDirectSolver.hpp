/*
 * Copyright (c) 2023 Nils Kohl.
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

#include "hyteg/eigen/EigenSparseMatrix.hpp"
#include "hyteg/eigen/EigenVector.hpp"
#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/sparseassembly/DirichletBCs.hpp"

namespace hyteg {

/// \brief Wrapper for the sparse LU solver provided by the Eigen library.
///
/// This solver should only be used for single-core testing purposes.
/// It is not possible to run this solver in an MPI-parallel fashion.
template < class OperatorType >
class EigenSparseDirectSolver : public Solver< OperatorType >
{
 public:
   using FunctionType = typename OperatorType::srcType;

   EigenSparseDirectSolver( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t& level )
   : storage_( storage )
   , level_( level )
   , numerator_( typename OperatorType::srcType::template FunctionType< idx_t >( "numerator", storage, level, level ) )
   , hasBeenAssembled_( false )
   , reassembleMatrix_( false )
   {
      numerator_.enumerate( level );
   }

   /// \brief If set to true, the operator is reassembled for every solve / manual assembly call.
   ///        Default is false.
   void setReassembleMatrix( bool reassembleMatrix ) { reassembleMatrix_ = reassembleMatrix; }

   /// \brief Assembles and factorizes the sparse matrix if either the matrix has not been assembled previously
   /// or reassembly is forced.
   void assembleAndFactorize( const OperatorType& A )
   {
      if ( hasBeenAssembled_ && !reassembleMatrix_ )
      {
         return;
      }

      storage_->getTimingTree()->start( "Matrix assembly" );

      mat_ = createEigenSparseMatrixFromOperator( A, level_, numerator_, numerator_ );

      std::vector< idx_t > bcIndices;
      applyDirichletBC( numerator_, bcIndices, level_ );

      for ( auto row : bcIndices )
      {
         for ( Eigen::SparseMatrix< real_t, Eigen::RowMajor >::InnerIterator it( mat_, static_cast< Eigen::Index >( row ) ); it;
               ++it )
         {
            if ( it.col() == row )
            {
               it.valueRef() = 1.0;
            }
            else
            {
               it.valueRef() = 0.0;
            }
         }
      }

      storage_->getTimingTree()->stop( "Matrix assembly" );

      storage_->getTimingTree()->start( "Factorization" );

      solver_.compute( mat_ );
      WALBERLA_CHECK_EQUAL( solver_.info(), Eigen::Success, "Eigen sparse LU factorization failed." )

      storage_->getTimingTree()->stop( "Factorization" );
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level )
   {
      WALBERLA_CHECK_EQUAL( level, level_ );

      walberla::WcTimer timer;

      storage_->getTimingTree()->start( "Eigen LU Solver" );
      storage_->getTimingTree()->start( "Setup" );
      storage_->getTimingTree()->start( "matrix setup" );

      // If the numerator was constructed internally we should copy the Boundary Condition info to it, because that is the only
      // way to get it into hyteg::applyDirichletBC() where the corresponding DoF indices will be computed
      numerator_.copyBoundaryConditionFromFunction( x );

      assembleAndFactorize( A );

      storage_->getTimingTree()->stop( "matrix setup" );

      storage_->getTimingTree()->start( "RHS vector setup" );

      rhs_ = createEigenVectorFromFunction( b, numerator_, level );
      assignEigenVectorFromFunction( x, rhs_, numerator_, level, DirichletBoundary );

      storage_->getTimingTree()->stop( "RHS vector setup" );

      storage_->getTimingTree()->stop( "Setup" );

      storage_->getTimingTree()->start( "Solver" );

      VectorXr sol = solver_.solve( rhs_ );
      WALBERLA_CHECK_EQUAL( solver_.info(), Eigen::Success, "Eigen sparse LU solve failed." )

      storage_->getTimingTree()->stop( "Solver" );

      assignFunctionFromEigenVector( sol, x, numerator_, level );

      storage_->getTimingTree()->stop( "Eigen LU Solver" );
   }

 private:
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              level_;

   Eigen::SparseLU< Eigen::SparseMatrix< real_t > > solver_;

   Eigen::SparseMatrix< real_t, Eigen::RowMajor >                 mat_;
   VectorXr                                                       rhs_;
   typename OperatorType::srcType::template FunctionType< idx_t > numerator_;

   bool hasBeenAssembled_;
   bool reassembleMatrix_;
};

} // namespace hyteg
