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
#include "hyteg/geometry/Polygons.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::int_c;
using walberla::real_c;
using walberla::math::pi;

namespace hyteg {

class TokamakMap : public GeometryMap
{
 public:
   TokamakMap( const Cell&                  cell,
               const SetupPrimitiveStorage& setupStorage,
               uint_t                       numSlices,
               uint_t                       numRadialEdges,
               real_t                       innerRadius,
               real_t                       outerRadius,
               real_t                       radiusZ,
               uint_t                       cutSide,
               uint_t                       cutTopAndBottom,
               real_t                       blendingCenterRadius )
   : numSlices_( numSlices )
   , numRadialEdges_( numRadialEdges )
   , innerRadius_( innerRadius )
   , outerRadius_( outerRadius )
   , radiusZ_( radiusZ )
   , cutSide_( cutSide )
   , cutTopAndBottom_( cutTopAndBottom )
   , blendingCenterRadius_( blendingCenterRadius )
   {
      identifyPrism( cell );
   }

   TokamakMap( const Face&                  face,
               const SetupPrimitiveStorage& setupStorage,
               uint_t                       numSlices,
               uint_t                       numRadialEdges,
               real_t                       innerRadius,
               real_t                       outerRadius,
               real_t                       radiusZ,
               uint_t                       cutSide,
               uint_t                       cutTopAndBottom,
               real_t                       blendingCenterRadius )
   : numSlices_( numSlices )
   , numRadialEdges_( numRadialEdges )
   , innerRadius_( innerRadius )
   , outerRadius_( outerRadius )
   , radiusZ_( radiusZ )
   , cutSide_( cutSide )
   , cutTopAndBottom_( cutTopAndBottom )
   , blendingCenterRadius_( blendingCenterRadius )
   {
      std::vector< PrimitiveID > neighborCells;
      face.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TokamakMap( const Edge&                  edge,
               const SetupPrimitiveStorage& setupStorage,
               uint_t                       numSlices,
               uint_t                       numRadialEdges,
               real_t                       innerRadius,
               real_t                       outerRadius,
               real_t                       radiusZ,
               uint_t                       cutSide,
               uint_t                       cutTopAndBottom,
               real_t                       blendingCenterRadius )
   : numSlices_( numSlices )
   , numRadialEdges_( numRadialEdges )
   , innerRadius_( innerRadius )
   , outerRadius_( outerRadius )
   , radiusZ_( radiusZ )
   , cutSide_( cutSide )
   , cutTopAndBottom_( cutTopAndBottom )
   , blendingCenterRadius_( blendingCenterRadius )
   {
      std::vector< PrimitiveID > neighborCells;
      edge.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TokamakMap( const Vertex&                vertex,
               const SetupPrimitiveStorage& setupStorage,
               uint_t                       numSlices,
               uint_t                       numRadialEdges,
               real_t                       innerRadius,
               real_t                       outerRadius,
               real_t                       radiusZ,
               uint_t                       cutSide,
               uint_t                       cutTopAndBottom,
               real_t                       blendingCenterRadius )
   : numSlices_( numSlices )
   , numRadialEdges_( numRadialEdges )
   , innerRadius_( innerRadius )
   , outerRadius_( outerRadius )
   , radiusZ_( radiusZ )
   , cutSide_( cutSide )
   , cutTopAndBottom_( cutTopAndBottom )
   , blendingCenterRadius_( blendingCenterRadius )
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
      auto innerSliceNormalD    = -innerPrismMidPoint_.dot( innerPrismNormal_ );
      auto distanceToInnerPlane = innerPrismNormal_.dot( xold ) + innerSliceNormalD;
      auto prismRadialThickness = ( outerPrismMidPoint_ - innerPrismMidPoint_ ).norm();

      auto newRadialDistance =
          innerRadius_ + ( distanceToInnerPlane / prismRadialThickness ) * ( outerRadiusRespectingCut_ - innerRadius_ );
      auto phi = std::atan2( xold[1], xold[0] );
      phi      = std::fmod( ( phi + 2 * pi ), 2 * pi );

      xnew[0] = newRadialDistance * std::cos( phi );
      xnew[1] = newRadialDistance * std::sin( phi );
      xnew[2] = xold[2];

      // identify all vertices of the slice-polygon
      auto numTangentialEdgesRespectingCut = numRadialEdges_ - cutTopAndBottom_;
      auto radialEdgeLength                = ( outerRadius_ - innerRadius_ ) / real_c( numRadialEdges_ );
      auto tangentialEdgeLength            = radiusZ_ / real_c( numRadialEdges_ );

      Point3D polygonVertexCenterOuter(
          { outerRadiusRespectingCut_ * std::cos( phi ), outerRadiusRespectingCut_ * std::sin( phi ), 0 } );
      Point3D polygonVertexCenterCutTop =
          polygonVertexCenterOuter + Point3D( { 0, 0, real_c( cutSide_ ) * tangentialEdgeLength } );
      Point3D polygonVertexCenterCutBot =
          polygonVertexCenterOuter + Point3D( { 0, 0, -real_c( cutSide_ ) * tangentialEdgeLength } );

      Point3D polygonVertexTop( { innerRadius_ * std::cos( phi ),
                                  innerRadius_ * std::sin( phi ),
                                  real_c( numTangentialEdgesRespectingCut ) * tangentialEdgeLength } );
      Point3D polygonVertexTopCut = polygonVertexTop + Point3D( { real_c( cutTopAndBottom_ ) * radialEdgeLength * std::cos( phi ),
                                                                  real_c( cutTopAndBottom_ ) * radialEdgeLength * std::sin( phi ),
                                                                  0 } );

      Point3D polygonVertexBot( { innerRadius_ * std::cos( phi ),
                                  innerRadius_ * std::sin( phi ),
                                  -real_c( numTangentialEdgesRespectingCut ) * tangentialEdgeLength } );
      Point3D polygonVertexBotCut = polygonVertexBot + Point3D( { real_c( cutTopAndBottom_ ) * radialEdgeLength * std::cos( phi ),
                                                                  real_c( cutTopAndBottom_ ) * radialEdgeLength * std::sin( phi ),
                                                                  0 } );

      std::vector< Point3D > slicePolygon;

      slicePolygon.push_back( polygonVertexCenterOuter );
      if ( cutSide_ )
      {
         slicePolygon.push_back( polygonVertexCenterCutTop );
      }
      if ( cutTopAndBottom_ )
      {
         slicePolygon.push_back( polygonVertexTopCut );
      }
      slicePolygon.push_back( polygonVertexTop );
      slicePolygon.push_back( polygonVertexBot );

      if ( cutTopAndBottom_ )
      {
         slicePolygon.push_back( polygonVertexBotCut );
      }
      if ( cutSide_ )
      {
         slicePolygon.push_back( polygonVertexCenterCutBot );
      }

      auto center = blendingCenterRadius_ * Point3D( { std::cos( phi ), std::sin( phi ), 0 } );

      real_t r;
      real_t theta;

      fractionalRadiusToPolygonBoundary( xnew, center, slicePolygon, r, theta );

      //      auto tubeRadius  = 0.5;
      //      auto torusRadius = 0.8;

      //      xnew[0] = ( torusRadius + r * tubeRadius * std::cos( angle ) ) * std::cos( phi );
      //      xnew[1] = ( torusRadius + r * tubeRadius * std::cos( angle ) ) * std::sin( phi );
      //      xnew[2] = r * tubeRadius * std::sin( angle );

      auto r_z  = r * std::sin( theta );
      auto r_xy = r * std::cos( theta );

      real_t r0           = outerRadius_;
      real_t r1           = innerRadius_;
      real_t r2           = r1 * 1.55;
      real_t arcsin_delta = 0.5;

      xnew[0] = ( 1 + r * ( r1 / r0 ) * std::cos( theta + arcsin_delta * std::sin( theta ) ) ) * std::cos( phi );
      xnew[1] = ( 1 + r * ( r1 / r0 ) * std::cos( theta + arcsin_delta * std::sin( theta ) ) ) * std::sin( phi );
      xnew[2] = r * ( r2 / r0 ) * std::sin( theta );
   }

   real_t evalDF( const Point3D& x, Matrix3r& DFx ) const final
   {
      WALBERLA_ABORT( "Not implemented." )
      return 0;
   }

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
                       uint_t                 cutSide,
                       uint_t                 cutTopAndBottom,
                       real_t                 blendingCenterRadius )
   {
      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(),
                                      std::make_shared< TokamakMap >( cell,
                                                                      setupStorage,
                                                                      numSlices,
                                                                      numRadialEdges,
                                                                      innerRadius,
                                                                      outerRadius,
                                                                      radiusZ,
                                                                      cutSide,
                                                                      cutTopAndBottom,
                                                                      blendingCenterRadius ) );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(),
                                      std::make_shared< TokamakMap >( face,
                                                                      setupStorage,
                                                                      numSlices,
                                                                      numRadialEdges,
                                                                      innerRadius,
                                                                      outerRadius,
                                                                      radiusZ,
                                                                      cutSide,
                                                                      cutTopAndBottom,
                                                                      blendingCenterRadius ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(),
                                      std::make_shared< TokamakMap >( edge,
                                                                      setupStorage,
                                                                      numSlices,
                                                                      numRadialEdges,
                                                                      innerRadius,
                                                                      outerRadius,
                                                                      radiusZ,
                                                                      cutSide,
                                                                      cutTopAndBottom,
                                                                      blendingCenterRadius ) );
      }

      for ( auto it : setupStorage.getVertices() )
      {
         Vertex& vertex = *it.second;
         setupStorage.setGeometryMap( vertex.getID(),
                                      std::make_shared< TokamakMap >( vertex,
                                                                      setupStorage,
                                                                      numSlices,
                                                                      numRadialEdges,
                                                                      innerRadius,
                                                                      outerRadius,
                                                                      radiusZ,
                                                                      cutSide,
                                                                      cutTopAndBottom,
                                                                      blendingCenterRadius ) );
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
   void identifyPrism( const Cell& cell )
   {
      auto    coords = cell.getCoordinates();
      Point3D centroid( { 0, 0, 0 } );
      for ( uint_t i = 0; i < 4; i++ )
      {
         centroid += coords[i];
      }
      centroid *= 0.25;

      auto radialEdgeLength     = ( outerRadius_ - innerRadius_ ) / real_c( numRadialEdges_ );
      auto tangentialEdgeLength = radiusZ_ / real_c( numRadialEdges_ );

      auto phi = std::atan2( centroid[1], centroid[0] );
      phi      = std::fmod( ( phi + 2 * pi ), 2 * pi );

      auto prismAngle = ( 2 * pi ) / real_c( numSlices_ );
      auto prismID_   = uint_c( phi / prismAngle );

      auto leftSlicePhi  = real_c( prismID_ ) * prismAngle;
      auto rightSlicePhi = real_c( ( prismID_ + 1 ) % numSlices_ ) * prismAngle;

      Point3D innerSliceLeftPoint( { innerRadius_ * std::cos( leftSlicePhi ), innerRadius_ * std::sin( leftSlicePhi ), 0 } );
      Point3D innerSliceRightPoint( { innerRadius_ * std::cos( rightSlicePhi ), innerRadius_ * std::sin( rightSlicePhi ), 0 } );
      innerPrismMidPoint_ = innerSliceLeftPoint + 0.5 * ( innerSliceRightPoint - innerSliceLeftPoint );
      innerPrismNormal_   = innerPrismMidPoint_ / innerPrismMidPoint_.norm();

      outerRadiusRespectingCut_ = outerRadius_ - real_c( cutSide_ ) * radialEdgeLength;

      Point3D outerEdgeLeftPoint(
          { outerRadiusRespectingCut_ * std::cos( leftSlicePhi ), outerRadiusRespectingCut_ * std::sin( leftSlicePhi ), 0 } );
      Point3D outerEdgeRightPoint(
          { outerRadiusRespectingCut_ * std::cos( rightSlicePhi ), outerRadiusRespectingCut_ * std::sin( rightSlicePhi ), 0 } );
      outerPrismMidPoint_ = outerEdgeLeftPoint + 0.5 * ( outerEdgeRightPoint - outerEdgeLeftPoint );

      //      ///////////////
      //      /// Sectors ///
      //      ///////////////
      //
      //      // define a center, could be chosen "arbitrarily"
      //      // blendingCenterRadius_ = innerRadius_ + 0.5 * ( outerRadiusRespectingCut_ - innerRadius_ );
      //
      //      Point2D center2D( { blendingCenterRadius_, 0 } );
      //
      //      Point2D bottomOuter = Point2D( { outerRadiusRespectingCut_, 0 } );
      //      Point2D outerCut    = bottomOuter + Point2D( { 0, tangentialEdgeLength } );
      //
      //      Point2D topCut( { innerRadius_ + radialEdgeLength, real_c( numRadialEdges_ - 1 ) * radialEdgeLength } );
      //      Point2D top = topCut + Point2D( { -radialEdgeLength, tangentialEdgeLength } );
      //      if ( cutTopAndBottom_ )
      //      {
      //         top -= Point2D( { 0, tangentialEdgeLength } );
      //      }
      //
      //      Point2D innerBottom( { innerRadius_, 0 } );
      //
      //      bottomOuter -= center2D;
      //      outerCut -= center2D;
      //      topCut -= center2D;
      //      top -= center2D;
      //      innerBottom -= center2D;
      //
      //      // angle = atan2(vector2.y, vector2.x) - atan2(vector1.y, vector1.x);
      //
      //      sectorOuterRadiiFromCenter_.push_back( bottomOuter.norm() );
      //      sectorOuterAnglesFromCenter_.push_back( std::atan2( bottomOuter[1], bottomOuter[0] ) );
      //
      //      if ( cutSide_ )
      //      {
      //         auto v1 = -bottomOuter;
      //         auto v2 = outerCut - bottomOuter;
      //         sectorOuterOpposingAnglesFromCenter_.push_back( std::atan2( v2[1], v2[0] ) - std::atan2( v1[1], v1[0] ) );
      //      }
      //      else
      //      {
      //         auto v1 = -bottomOuter;
      //         auto v2 = top - bottomOuter;
      //         sectorOuterOpposingAnglesFromCenter_.push_back( std::atan2( v2[1], v2[0] ) - std::atan2( v1[1], v1[0] ) );
      //      }
      //
      //      if ( cutSide_ )
      //      {
      //         sectorOuterRadiiFromCenter_.push_back( outerCut.norm() );
      //         sectorOuterAnglesFromCenter_.push_back( std::atan2( outerCut[1], outerCut[0] ) );
      //         auto v1 = -outerCut;
      //         auto v2 = topCut - outerCut;
      //         sectorOuterOpposingAnglesFromCenter_.push_back( std::atan2( v2[1], v2[0] ) - std::atan2( v1[1], v1[0] ) );
      //
      //         sectorOuterRadiiFromCenter_.push_back( topCut.norm() );
      //         sectorOuterAnglesFromCenter_.push_back( std::atan2( topCut[1], topCut[0] ) );
      //         v1 = -topCut;
      //         v2 = top - topCut;
      //         sectorOuterOpposingAnglesFromCenter_.push_back( std::atan2( v2[1], v2[0] ) - std::atan2( v1[1], v1[0] ) );
      //      }
      //
      //      sectorOuterRadiiFromCenter_.push_back( top.norm() );
      //      sectorOuterAnglesFromCenter_.push_back( std::atan2( top[1], top[0] ) );
      //      auto v1 = -top;
      //      auto v2 = innerBottom - top;
      //      sectorOuterOpposingAnglesFromCenter_.push_back( std::atan2( v2[1], v2[0] ) - std::atan2( v1[1], v1[0] ) );
      //
      //      sectorOuterRadiiFromCenter_.push_back( innerBottom.norm() );
   }

   // number of slices or prisms
   uint_t numSlices_;

   // number of edges in radial direction, including the tip, even if cutSide == true
   uint_t numRadialEdges_;

   // inner and outer radii at the slices that interface two prisms, also count triangle tip if cutSide == true
   real_t innerRadius_;
   real_t outerRadius_;

   // height in z-direction from z == 0, including tip even if cutTopAndBottom == true
   real_t radiusZ_;

   uint_t cutSide_;
   uint_t cutTopAndBottom_;

   // midpoint of the inner prism plane, before blending!
   Point3D innerPrismMidPoint_;
   // outward pointing normal at the inner prism mid point
   Point3D innerPrismNormal_;
   // midpoint outer prism line (or plane if cutSide == true), before blending!
   Point3D outerPrismMidPoint_;

   // outer radius of the slice taking the cutSide setting into account, if cutSide == false, this is equal to outerRadius_ of the slice
   real_t outerRadiusRespectingCut_;

   // radius of the center point to perform the blending of the slices, assuming the center point is always at z == 0
   real_t blendingCenterRadius_;
   //
   //   std::vector< real_t > sectorOuterRadiiFromCenter_;
   //   std::vector< real_t > sectorOuterAnglesFromCenter_;
   //   std::vector< real_t > sectorOuterOpposingAnglesFromCenter_;
};

} // end of namespace hyteg
