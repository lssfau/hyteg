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

#include "core/DataTypes.h"

#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/solvers/Solver.hpp"

namespace hyteg {

template < class OperatorType >
class WeightedJacobiSmoother : public Solver< OperatorType >
{
 public:
   using FunctionType = typename OperatorType::srcType;
   using ValueType    = typename FunctionTrait< FunctionType >::ValueType;

   WeightedJacobiSmoother( const std::shared_ptr< PrimitiveStorage >& storage,
                           uint_t                                     minLevel,
                           uint_t                                     maxLevel,
                           const real_t&                              relax )
   : relax_( relax )
   , tmp_( "tmp_weighted_jacobi", storage, minLevel, maxLevel )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary )
   {}

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      x.getStorage()->getTimingTree()->start( "Weighted Jacobi" );
      tmp_.copyBoundaryConditionFromFunction( x );
      tmp_.assign( { walberla::numeric_cast< ValueType >( 1.0 ) }, { x }, level, All );
      if ( const auto* A_jac = dynamic_cast< const WeightedJacobiSmoothable< typename OperatorType::srcType >* >( &A ) )
      {
         A_jac->smooth_jac( x, b, tmp_, relax_, level, flag_ );
      }
      else
      {
         throw std::runtime_error( "The WeightedJacobiSmoother requires the WeightedJacobiSmoothable interface." );
      }
      x.getStorage()->getTimingTree()->stop( "Weighted Jacobi" );
   }

 private:
   real_t                         relax_;
   typename OperatorType::srcType tmp_;
   DoFType                        flag_;
};

} // namespace hyteg
