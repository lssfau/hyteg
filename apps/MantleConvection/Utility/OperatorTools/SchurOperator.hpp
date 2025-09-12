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

#include "hyteg/operators/Operator.hpp"

namespace MantleConvection {

template < class PressureFunctionType >
class SchurOperator : public hyteg::Operator< PressureFunctionType, PressureFunctionType >
{
 public:
   SchurOperator( const std::shared_ptr< hyteg::PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : hyteg::Operator< PressureFunctionType, PressureFunctionType >( storage, minLevel, maxLevel )
   {}

   void apply( const PressureFunctionType& src,
               const PressureFunctionType& dst,
               const uint_t                level,
               const hyteg::DoFType        flag,
               const hyteg::UpdateType     updateType = hyteg::Replace ) const
   {
      WALBERLA_UNUSED( src );
      WALBERLA_UNUSED( dst );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_UNUSED( updateType );

      WALBERLA_ABORT(
          "The SchurOperator class is not meant to be used but as a stand-in for the operator type that Schur complement solvers are solving for." )
   }
};

} // namespace MantleConvection
