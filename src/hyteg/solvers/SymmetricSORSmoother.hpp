/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {

template < class OperatorType >
class SymmetricSORSmoother : public Solver< OperatorType >
{
 public:
   SymmetricSORSmoother( const real_t& relax )
   : relax_( relax )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary )
   {}

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      if ( const auto* A_sor = dynamic_cast< const SORSmoothable< typename OperatorType::srcType >* >( &A ) )
      {
         A_sor->smooth_sor( x, b, relax_, level, flag_ );
      }
      else
      {
         throw std::runtime_error( "The symmetric SOR-Operator requires the SORSmoothable interface." );
      }

      if ( const auto* A_sor_bw = dynamic_cast< const SORBackwardsSmoothable< typename OperatorType::srcType >* >( &A ) )
      {
         A_sor_bw->smooth_sor_backwards( x, b, relax_, level, flag_ );
      }
      else
      {
         throw std::runtime_error( "The symmetric SOR-Operator requires the SORBackwardsSmoothable interface." );
      }
   }

 private:
   real_t  relax_;
   DoFType flag_;
};

} // namespace hyteg
