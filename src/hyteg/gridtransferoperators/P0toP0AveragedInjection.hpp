/*
 * Copyright (c) 2025 Ponsuganth Ilangovan.
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

#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/types/Averaging.hpp"

namespace hyteg {

class P0toP0AveragedInjection : public RestrictionOperator< P0Function< real_t > >
{
 public:
   /// \brief Transfers a P0 function from sourceLevel to lower level(s)
   ///
   /// The transfer happens with a user-specified averaging method from hyteg::AveragingType
   /// Note that, some averaging methods from hyteg::AveragingType are not implemented
   /// In 3D (2D), values from the 8 (4) child tetrahedrons (triangles) are averaged
   /// and is specified to their parent tetrahedron (triangle)
   /// In addition if a volume weighted averaging has to be done 
   /// can also be specified
   ///
   /// \param averagingType  user-specified averaging method from hyteg::AveragingType
   /// \param volumeWeighted if the averaging needs to be volume weighted
   ///
   P0toP0AveragedInjection( hyteg::AveragingType averagingType, bool volumeWeighted )
   : averagingType_( averagingType )
   , volumeWeighted_( volumeWeighted )
   {}

   /// \brief Transfers a P0 function from sourceLevel to the next lower level
   ///
   /// \param function    The P0 function that needs to be average-transferred
   /// \param sourceLevel Level from which the function should be considered for averaging
   /// \param flag        This parameter is not used at the moment
   ///
   void restrict( const P0Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag = All ) const override;
   
   /// \brief Transfers a P0 function from sourceLevel to all lower levels 
   ///        till the function's minLevel
   ///
   /// \param function    The P0 function that needs to be average-transferred
   /// \param sourceLevel Level from which the function should be considered for averaging
   ///
   void restrictToAllLowerLevels( const P0Function< real_t >& function, const uint_t& sourceLevel );

 private:
   hyteg::AveragingType averagingType_;
   bool                 volumeWeighted_;
};

} // namespace hyteg
