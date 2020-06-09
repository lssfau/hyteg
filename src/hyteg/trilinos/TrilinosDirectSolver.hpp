/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "hyteg/solvers/Solver.hpp"
#include "hyteg/trilinos/Amesos2Wrapper.hpp"
#include "hyteg/trilinos/TrilinosSparseMatrix.hpp"
#include "hyteg/trilinos/TrilinosVector.hpp"

namespace hyteg {
namespace trilinos {

using walberla::uint_t;

using Teuchos::RCP;
using Teuchos::rcp;

enum class TrilinosDirectSolverType
{
   MUMPS
};

/// \brief Solver interface that wraps the Trilinos package Amesos2 that implements various sparse direct solvers.
template < typename OperatorType >
class TrilinosDirectSolver : public Solver< OperatorType >
{
 public:
   template < typename ValueType >
   using FunctionTemplate = typename OperatorType::srcType::template FunctionType< ValueType >;

   TrilinosDirectSolver( const TrilinosDirectSolverType&            solverType,
                         const std::shared_ptr< PrimitiveStorage >& storage,
                         const uint_t&                              level,
                         const DoFType&                             flag )
   : solverType_( solverType )
   , level_( level )
   , ATrilinos_( storage, level )
   , xTrilinos_( storage, level )
   , bTrilinos_( storage, level )
   , numerator_( "trilinos_numerator", storage, level, level )
   , flag_( flag )
   {
      numerator_.enumerate( level );
      std::string solverTypeString = "MUMPS";

      solver_ = Amesos2::create< typename TrilinosSparseMatrix< OperatorType, FunctionTemplate >::MatrixType,
                                 typename TrilinosVector< FunctionTemplate >::VectorType >
                ( solverTypeString, ATrilinos_.getTpetraMatrix(), xTrilinos_.getTpetraVector(), bTrilinos_.getTpetraVector() );
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const uint_t                          level )
   {
      WALBERLA_CHECK_EQUAL( level, level_, "Cannot employ Trilinos direct solver on level " << level );

      // assemble the sparse system
      ATrilinos_.assembleSystem( A, numerator_, flag_ );
      ATrilinos_.applyDirichletBoundaryConditions( numerator_ );

      // copy the vectors
      bTrilinos_.fillFromFunction( b, numerator_, flag_ );
      bTrilinos_.fillFromFunction( x, numerator_, DirichletBoundary );
      xTrilinos_.fillFromFunction( x, numerator_ );

      solver_->symbolicFactorization().numericFactorization().solve();

      // write back
      xTrilinos_.writeToFunction( x, numerator_ );
   }

 private:
   TrilinosDirectSolverType solverType_;
   uint_t                   level_;

   TrilinosSparseMatrix< OperatorType, FunctionTemplate > ATrilinos_;
   TrilinosVector< FunctionTemplate >                     xTrilinos_;
   TrilinosVector< FunctionTemplate >                     bTrilinos_;
   FunctionTemplate< PetscInt >                           numerator_;
   DoFType                                                flag_;
   RCP< Amesos2::Solver< typename TrilinosSparseMatrix< OperatorType, FunctionTemplate >::MatrixType,
                         typename TrilinosVector< FunctionTemplate >::VectorType > >
       solver_;
};

} // namespace trilinos
} // namespace hyteg