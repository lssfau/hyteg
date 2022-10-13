/*
 * Copyright (c) 2022 Daniel Bauer.
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
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"

namespace hyteg {
namespace n1e1 {

class N1E1toN1E1Restriction : public RestrictionOperator< N1E1VectorFunction< real_t > >
{
 public:
   void restrict( const N1E1VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      restrictAdditively( function, sourceLevel, flag, Replace );
   }

 private:
   void restrictAdditively( const N1E1VectorFunction< real_t >& function,
                            const uint_t&                       sourceLevel,
                            const DoFType&                      flag,
                            const UpdateType&                   updateType ) const;

   void setToZeroIfNeededMacroEdge( FunctionMemory< real_t >* data, const uint_t& level, const UpdateType& updateType ) const;
   void setToZeroIfNeededMacroFace( FunctionMemory< real_t >* data, const uint_t& level, const UpdateType& updateType ) const;
   void setToZeroIfNeededMacroCell( FunctionMemory< real_t >* data, const uint_t& level, const UpdateType& updateType ) const;

   void restrictMacroEdge( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const;
   void restrictMacroFace( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const;
   void restrictMacroCell( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const;
};

} // namespace n1e1
} // namespace hyteg
