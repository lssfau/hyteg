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

#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"

namespace hyteg {
namespace n1e1 {

class N1E1toN1E1Prolongation : public ProlongationOperator< N1E1VectorFunction< real_t > >
{
 public:
   void prolongate( const N1E1VectorFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      prolongateAdditively( function, sourceLevel, flag, Replace );
   }
   void prolongateAndAdd( const N1E1VectorFunction< real_t >& function,
                          const walberla::uint_t&             sourceLevel,
                          const DoFType&                      flag ) const override
   {
      prolongateAdditively( function, sourceLevel, flag, Add );
   }

 private:
   void prolongateAdditively( const N1E1VectorFunction< real_t >& function,
                              const uint_t&                       sourceLevel,
                              const DoFType&                      flag,
                              const UpdateType&                   updateType ) const;

   void prolongateMacroEdge( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const;
   void prolongateMacroFace( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const;
   void prolongateMacroCell( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const;
};

} // namespace n1e1
} // namespace hyteg
