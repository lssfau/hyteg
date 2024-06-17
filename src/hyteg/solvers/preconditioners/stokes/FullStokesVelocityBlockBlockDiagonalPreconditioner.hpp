/*
 * Copyright (c) 2017-2020 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/EmptySolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

namespace hyteg {

/// \brief Applies a Solver to the velocity block for a Stokes system.
///
/// This preconditioner can be used in the UzawaSmoother to operate on the velocity block
/// of the Stokes system.
///
template < class OperatorType >
class FullStokesVelocityBlockBlockDiagonalPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   FullStokesVelocityBlockBlockDiagonalPreconditioner(
       const std::shared_ptr< PrimitiveStorage >&                                    storage,
       const std::shared_ptr< Solver< typename OperatorType::VelocityOperator_T > >& scalarVelocityPreconditioner )
   : hasGlobalCells_( storage->hasGlobalCells() )
   , scalarVelocityPreconditioner_( scalarVelocityPreconditioner )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, uint_t level ) override
   {
      scalarVelocityPreconditioner_->solve( A.getA(), x.uvw(), b.uvw(), level );
   }

 private:
   bool                                                                   hasGlobalCells_;
   std::shared_ptr< Solver< typename OperatorType::VelocityOperator_T > > scalarVelocityPreconditioner_;
};

} // namespace hyteg
