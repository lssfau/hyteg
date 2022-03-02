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

#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/Solver.hpp"

namespace hyteg {

template < class OperatorType, class pressureBlockPreconditionerType >
class StokesPressureBlockPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   StokesPressureBlockPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : pressureBlockPreconditioner_( std::make_shared< pressureBlockPreconditionerType >( storage, minLevel, maxLevel ) )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   {}

   void solve( const OperatorType&, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      x.assign( {1.0}, {b}, level, flag_ );
      pressureBlockPreconditioner_->apply( b.p(), x.p(), level, flag_, Replace );
   }

 private:
   std::shared_ptr< pressureBlockPreconditionerType > pressureBlockPreconditioner_;
   hyteg::DoFType                                       flag_;
};

} // namespace hyteg
