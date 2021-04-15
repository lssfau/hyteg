/*
 * Copyright (c) 2017-2021 Nils Kohl.
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

#include <cmath>

#include "core/math/Constants.h"

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::int_c;
using walberla::real_c;
using walberla::math::pi;

namespace hyteg {

/// The class implements a generic affine mapping in 3D
///
/// The affine mapping is characterised by a matrix \f$M\f$ and a vector \f$v\f$
/// and defined as
/// \f[
/// x \mapsto M x + v
/// \f]
class TokamakMap : public GeometryMap
{
 public:

   TokamakMap( const Cell& cell,
               const SetupPrimitiveStorage& setupStorage,
               uint_t numSlices,
               uint_t numRadialEdges,
               real_t innerRadius,
               real_t outerRadius,
               real_t radiusZ,
               bool   cutSide,
               bool   cutTopAndBottom )
   : numSlices_( numSlices )
   , numRadialEdges_( numRadialEdges )
   , innerRadius_( innerRadius )
   , outerRadius_( outerRadius )
   , radiusZ_( radiusZ )
   , cutSide_( cutSide )
   , cutTopAndBottom_( cutTopAndBottom )
   {
      identifyPrism( cell );
   }

   TokamakMap( const Face& face,
               const SetupPrimitiveStorage& setupStorage,
               uint_t numSlices,
               uint_t numRadialEdges,
               real_t innerRadius,
               real_t outerRadius,
               real_t radiusZ,
               bool   cutSide,
               bool   cutTopAndBottom )
       : numSlices_( numSlices )
       , numRadialEdges_( numRadialEdges )
       , innerRadius_( innerRadius )
       , outerRadius_( outerRadius )
       , radiusZ_( radiusZ )
       , cutSide_( cutSide )
       , cutTopAndBottom_( cutTopAndBottom )
   {
      std::vector< PrimitiveID > neighborCells;
      face.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TokamakMap( const Edge& edge,
               const SetupPrimitiveStorage& setupStorage,
               uint_t numSlices,
               uint_t numRadialEdges,
               real_t innerRadius,
               real_t outerRadius,
               real_t radiusZ,
               bool   cutSide,
               bool   cutTopAndBottom )
       : numSlices_( numSlices )
       , numRadialEdges_( numRadialEdges )
       , innerRadius_( innerRadius )
       , outerRadius_( outerRadius )
       , radiusZ_( radiusZ )
       , cutSide_( cutSide )
       , cutTopAndBottom_( cutTopAndBottom )
   {
      std::vector< PrimitiveID > neighborCells;
      edge.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TokamakMap( const Vertex& vertex,
               const SetupPrimitiveStorage& setupStorage,
               uint_t numSlices,
               uint_t numRadialEdges,
               real_t innerRadius,
               real_t outerRadius,
               real_t radiusZ,
               bool   cutSide,
               bool   cutTopAndBottom )
       : numSlices_( numSlices )
       , numRadialEdges_( numRadialEdges )
       , innerRadius_( innerRadius )
       , outerRadius_( outerRadius )
       , radiusZ_( radiusZ )
       , cutSide_( cutSide )
       , cutTopAndBottom_( cutTopAndBottom )
   {
      std::vector< PrimitiveID > neighborCells;
      vertex.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TokamakMap( walberla::mpi::RecvBuffer& recvBuffer ) { WALBERLA_ABORT( "Deserialization not implemented for TokamakMap" ); }

   void evalF( const Point3D& xold, Point3D& xnew ) const
   {
      auto innerSliceNormalD = - innerSliceMidPoint_.dot( innerSliceNormal_ );
      auto distanceToInnerPlane = innerSliceNormal_.dot( xold ) + innerSliceNormalD;
      auto prismRadialThickness = (outerSliceMidPoint_ - innerSliceMidPoint_).norm();

      auto newRadialDistance = innerRadius_ + ( distanceToInnerPlane / prismRadialThickness ) * ( outerRadius_ - innerRadius_ );
      auto phi = std::atan2( xold[1], xold[0] );

      xnew[0] = newRadialDistance * std::cos( phi );
      xnew[1] = newRadialDistance * std::sin( phi );
      xnew[2] = xold[2];

      // all slices are equal now



   }

   real_t evalDF( const Point3D& x, Matrix3r& DFx ) const final { return 0; }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const
   {
      WALBERLA_ABORT( "Serialization not implemented for TokamakMap" );
   }

   static void setMap( SetupPrimitiveStorage& setupStorage,
                       uint_t                 numSlices,
                       uint_t                 numRadialEdges,
                       real_t                 innerRadius,
                       real_t                 outerRadius,
                       real_t                 radiusZ,
                       bool                   cutSide,
                       bool                   cutTopAndBottom )
   {
      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap(
             cell.getID(),
             std::make_shared< TokamakMap >(
                 cell, setupStorage, numSlices, numRadialEdges, innerRadius, outerRadius, radiusZ, cutSide, cutTopAndBottom ) );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap(
             face.getID(),
             std::make_shared< TokamakMap >(
                 face, setupStorage, numSlices, numRadialEdges, innerRadius, outerRadius, radiusZ, cutSide, cutTopAndBottom ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap(
             edge.getID(),
             std::make_shared< TokamakMap >(
                 edge, setupStorage, numSlices, numRadialEdges, innerRadius, outerRadius, radiusZ, cutSide, cutTopAndBottom ) );
      }

      for ( auto it : setupStorage.getVertices() )
      {
         Vertex& vertex = *it.second;
         setupStorage.setGeometryMap(
             vertex.getID(),
             std::make_shared< TokamakMap >(
                 vertex, setupStorage, numSlices, numRadialEdges, innerRadius, outerRadius, radiusZ, cutSide, cutTopAndBottom ) );
      }
   }

   /** @name 2D methods
   *    methods for 2D (class only provides a pseudo-implementation to satisfy requirements of base class)
   */
   ///@{
   void evalDF( const Point3D& x, Matrix2r& DFx ) const final
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFx );
      WALBERLA_ABORT( "TokamakMap::evalDF unimplemented for 2D!" );
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const final
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFinvx );
      WALBERLA_ABORT( "TokamakMap::evalDFinv unimplemented for 2D!" );
   }
   ///@}

 private:

   void identifyPrism( const Cell & cell )
   {
      auto coords = cell.getCoordinates();
      Point3D centroid( {0, 0, 0} );
      for ( uint_t i = 0; i < 4; i++ )
      {
         centroid += coords[i];
      }
      centroid *= 0.25;

      auto phi = std::atan2( centroid[1], centroid[0] ) + pi;
      auto prismAngle = (2 * pi) / real_c( numSlices_ );
      auto prismID_ = uint_c( phi / prismAngle );

      // find the three planes that enclose the prism with radial normals
      // anchor point can be set to anything for now
      auto leftSlicePhi = real_c( prismID_ ) * prismAngle;
      auto rightSlicePhi = real_c( (prismID_ + 1) % numSlices_ ) * prismAngle;

      Point3D innerSliceLeftPoint( { innerRadius_ * std::cos( leftSlicePhi ), innerRadius_ * std::sin( leftSlicePhi ), 0 } );
      Point3D innerSliceRightPoint( { innerRadius_ * std::cos( rightSlicePhi ), innerRadius_ * std::sin( rightSlicePhi ), 0 } );
      innerSliceMidPoint_ = innerSliceLeftPoint + 0.5 * ( innerSliceRightPoint - innerSliceLeftPoint );
      innerSliceNormal_ = innerSliceMidPoint_ / innerSliceMidPoint_.norm();

      Point3D outerEdgeLeftPoint( { outerRadius_ * std::cos( leftSlicePhi ), outerRadius_ * std::sin( leftSlicePhi ), 0 } );
      Point3D outerEdgeRightPoint( { outerRadius_ * std::cos( rightSlicePhi ), outerRadius_ * std::sin( rightSlicePhi ), 0 } );
      outerSliceMidPoint_ = outerEdgeLeftPoint + 0.5 * ( outerEdgeRightPoint - outerEdgeLeftPoint );

      WALBERLA_LOG_DEVEL_ON_ROOT( innerSliceMidPoint_ );

   }

   uint_t numSlices_;
   uint_t numRadialEdges_;
   real_t innerRadius_;
   real_t outerRadius_;
   real_t radiusZ_;
   bool   cutSide_;
   bool   cutTopAndBottom_;

   Point3D innerSliceMidPoint_;
   Point3D outerSliceMidPoint_;
   Point3D innerSliceNormal_;

};

} // end of namespace hyteg
