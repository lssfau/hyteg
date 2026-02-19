/*
 * Copyright (c) 2026 Marcus Mohr.
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

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::real_c;

template < typename ValueType >
class P3Function final : public Function< P3Function< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = P3Function< VType >;

   P3Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel );

   P3Function( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               BoundaryCondition                          boundaryCondition );

   uint_t getDimension() const override final { return 1u; }

   P3Function< ValueType >& operator[]( uint_t idx )
   {
      WALBERLA_ASSERT( idx == 0 );
      WALBERLA_UNUSED( idx );
      return *this;
   }

   const P3Function< ValueType >& operator[]( uint_t idx ) const
   {
      WALBERLA_ASSERT( idx == 0 );
      WALBERLA_UNUSED( idx );
      return *this;
   }

   inline const vertexdof::VertexDoFFunction< ValueType >&      getVertexDoFFunction() const { return vertexDoFFunction_; }
   inline const EdgeDoFFunction< ValueType >&                   getEdgeDoFFunctionRed() const { return edgeDoFFunctionRed_; }
   inline const EdgeDoFFunction< ValueType >&                   getEdgeDoFFunctionBlue() const { return edgeDoFFunctionBlue_; }
   inline const volumedofspace::VolumeDoFFunction< ValueType >& getFaceDoFFunction() const { return faceDoFFunction_; }

   /// Set all function DoFs to zero including the ones in the halos
   void setToZero( const uint_t level ) const override final
   {
      vertexDoFFunction_.setToZero( level );
      edgeDoFFunctionRed_.setToZero( level );
      edgeDoFFunctionBlue_.setToZero( level );
      faceDoFFunction_.setToZero( level );
   };

 private:
   using Function< P3Function< ValueType > >::communicators_;

   vertexdof::VertexDoFFunction< ValueType > vertexDoFFunction_;
   EdgeDoFFunction< ValueType >              edgeDoFFunctionBlue_;
   EdgeDoFFunction< ValueType >              edgeDoFFunctionRed_;

   // TODO: In 2D this is fine, but in 3D we need a true FaceDoFFunction!!!!
   volumedofspace::VolumeDoFFunction< ValueType > faceDoFFunction_;
};

extern template class P3Function< double >;
extern template class P3Function< int >;
extern template class P3Function< idx_t >;

} //namespace hyteg
