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
#include "core/config/Config.h"

#include "hyteg/functions/FunctionHistory.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

template < class DstFunctionType,
           class MassOperatorType,
           class AdditionalDataType = real_t,
           class SrcFunctionType    = DstFunctionType >
class TimeDiscretisationScheme
{
 public:
   virtual ~TimeDiscretisationScheme() = default;

   // Warning! Only apply this AFTER you have completely assembled the rest of the LHS!
   virtual void applyLHS( const FunctionHistory< SrcFunctionType, AdditionalDataType >& history,
                          const SrcFunctionType&                                        src,
                          const DstFunctionType&                                        dst,
                          const MassOperatorType&                                       massOperator,
                          const uint_t                                                  level,
                          const hyteg::DoFType                                          flag )
   {
      WALBERLA_UNUSED( history );
      WALBERLA_UNUSED( src );
      WALBERLA_UNUSED( dst );
      WALBERLA_UNUSED( massOperator );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );

      WALBERLA_ABORT( "applyLHS not implemented for TimeDiscretisationScheme abstract base class!" );
   }

   // Warning! Only apply this AFTER you have completely assembled the rest of the RHS!
   virtual void applyRHS( const FunctionHistory< SrcFunctionType, AdditionalDataType >& src,
                          const DstFunctionType&                                        dst,
                          const MassOperatorType&                                       massOperator,
                          const uint_t                                                  level,
                          const hyteg::DoFType                                          flag )
   {
      WALBERLA_UNUSED( src );
      WALBERLA_UNUSED( dst );
      WALBERLA_UNUSED( massOperator );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "applyRHS not implemented for TimeDiscretisationScheme abstract base class!" );
   }

   // Warning! Only apply this AFTER you have completely assembled the rest of the matrix!
   virtual void toMatrix( const FunctionHistory< SrcFunctionType, AdditionalDataType >&          src,
                          const MassOperatorType&                                                massOperator,
                          const std::shared_ptr< hyteg::SparseMatrixProxy >&                     mat,
                          const typename SrcFunctionType::template FunctionType< hyteg::idx_t >& srcIdx,
                          const typename DstFunctionType::template FunctionType< hyteg::idx_t >& dstIdx,
                          uint_t                                                                 level,
                          hyteg::DoFType                                                         flag )
   {
      WALBERLA_UNUSED( src );
      WALBERLA_UNUSED( massOperator );
      WALBERLA_UNUSED( mat );
      WALBERLA_UNUSED( srcIdx );
      WALBERLA_UNUSED( dstIdx );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "toMatrix not implemented for TimeDiscretisationScheme abstract base class!" );
   }
};

} // namespace hyteg