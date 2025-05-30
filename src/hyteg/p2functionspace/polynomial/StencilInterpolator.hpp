/*
* Copyright (c) 2017-2019 Benjamin Mann.
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

#include <array>
#include <hyteg/p2functionspace/polynomial/P2StencilPolynomial.hpp>
#include <hyteg/p2functionspace/variablestencil/P2VariableStencilCommon.hpp>
#include <hyteg/polynomial/LSQPInterpolator.hpp>

namespace hyteg {
namespace P2 {

template < LSQPType L, P2::NumStencilentries2D N >
class StencilInterpolator
{
 public:
   StencilInterpolator( uint_t polyDegree )
   {
      for ( uint_t i = 0; i < N; ++i )
      {
         data_[i] = std::make_shared< LSQPInterpolator< MonomialBasis2D, L > >( polyDegree );
      }
   }

   void interpolate( StencilPolynomial< N >& poly )
   {
      for ( uint_t i = 0; i < N; ++i )
      {
         data_[i]->interpolate( poly[i] );
      }
   }

   LSQPInterpolator< MonomialBasis2D, L >& operator[]( uint_t i ) { return *( data_[i] ); }

   const LSQPInterpolator< MonomialBasis2D, L >& operator[]( uint_t i ) const { return *( data_[i] ); }

 private:
   std::array< std::shared_ptr< LSQPInterpolator< MonomialBasis2D, L > >, N > data_;
};

} // namespace P2
} // namespace hyteg
