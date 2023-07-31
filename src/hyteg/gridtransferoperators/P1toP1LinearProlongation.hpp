/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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
#include "hyteg/p1functionspace/P1Function.hpp"

namespace hyteg {

template < typename ValueType = real_t >
class P1toP1LinearProlongation : public ProlongationOperator< P1Function< ValueType > >
{
 public:
   void prolongate( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      if ( function.getStorage()->hasGlobalCells() )
      {
         prolongate3DAdditively( function, sourceLevel, flag, Replace );
      }
      else
      {
         prolongate2DAdditively( function, sourceLevel, flag, Replace );
      }
   }

   void prolongateAndAdd( const P1Function< ValueType >& function,
                          const walberla::uint_t&        sourceLevel,
                          const DoFType&                 flag ) const override
   {
      if ( function.getStorage()->hasGlobalCells() )
      {
         prolongate3DAdditively( function, sourceLevel, flag, Add );
      }
      else
      {
         prolongate2DAdditively( function, sourceLevel, flag, Add );
      }
   }

 private:
   void prolongate2D( const P1Function< ValueType >& function, const uint_t& sourceLevel, const DoFType& flag ) const;
   void prolongate2DAdditively( const P1Function< ValueType >& function,
                                const uint_t&                  sourceLevel,
                                const DoFType&                 flag,
                                const UpdateType&              updateType ) const;

   void prolongate3DAdditively( const P1Function< ValueType >& function,
                                const uint_t&                  sourceLevel,
                                const DoFType&                 flag,
                                const UpdateType&              updateType ) const;

   void prolongateMacroVertex2D( const ValueType* src, ValueType* dst, const uint_t& sourceLevel ) const;

   void prolongateMacroEdge2D( const ValueType* src, ValueType* dst, const uint_t& sourceLevel ) const;

   void prolongateMacroFace2D( const ValueType* src, ValueType* dst, const uint_t& sourceLevel ) const;
};

} // namespace hyteg
