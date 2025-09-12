/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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
#include "core/Environment.h"
#include "core/config/Config.h"

#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/solvers/Solver.hpp"

namespace MantleConvection {

template < class OperatorType >
class ABlockSolver : public hyteg::Solver< OperatorType >
{
 public:
   virtual ~ABlockSolver() = default;

   ABlockSolver( std::string prefix = "" )
   : prefix_( prefix )
   {}

   virtual void solve( const OperatorType&                   A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const walberla::uint_t                level ) override
   {
      WALBERLA_UNUSED(A);
      WALBERLA_UNUSED(x);
      WALBERLA_UNUSED(b);
      WALBERLA_UNUSED(level);

      WALBERLA_ABORT( "solve not implemented for ABlockSolver abstract base class!" );
   };

   virtual std::ostream& print( std::ostream& os, uint_t offset = 0 ) const
   {
      WALBERLA_UNUSED( offset );
      WALBERLA_ABORT( "print not implemented for ABlockSolver abstract base class!" );
      return os;
   }

 protected:
   std::string prefix_;
};

} // namespace MantleConvection