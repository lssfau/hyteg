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
#include "hyteg/geometry/Torus.hpp"
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
   TokamakMap( const Cell&                  cell,
               const SetupPrimitiveStorage& setupStorage,
               uint_t                       numToroidalSlices,
               uint_t                       numPoloidalSlices,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r0,
               real_t                       r1,
               real_t                       r2 )
   : numToroidalSlices_( numToroidalSlices )
   , numPoloidalSlices_( numPoloidalSlices )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )
   , delta_( delta )
   , r0_( r0 )
   , r1_( r1 )
   , r2_( r2 )

   {
      identifyPrism( cell );
   }

   TokamakMap( const Face&                  face,
               const SetupPrimitiveStorage& setupStorage,
               uint_t                       numToroidalSlices,
               uint_t                       numPoloidalSlices,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r0,
               real_t                       r1,
               real_t                       r2 )
   : numToroidalSlices_( numToroidalSlices )
   , numPoloidalSlices_( numPoloidalSlices )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )
   , delta_( delta )
   , r0_( r0 )
   , r1_( r1 )
   , r2_( r2 )
   {
      std::vector< PrimitiveID > neighborCells;
      face.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TokamakMap( const Edge&                  edge,
               const SetupPrimitiveStorage& setupStorage,
               uint_t                       numToroidalSlices,
               uint_t                       numPoloidalSlices,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r0,
               real_t                       r1,
               real_t                       r2 )
   : numToroidalSlices_( numToroidalSlices )
   , numPoloidalSlices_( numPoloidalSlices )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )
   , delta_( delta )
   , r0_( r0 )
   , r1_( r1 )
   , r2_( r2 )
   {
      std::vector< PrimitiveID > neighborCells;
      edge.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TokamakMap( const Vertex&                vertex,
               const SetupPrimitiveStorage& setupStorage,
               uint_t                       numToroidalSlices,
               uint_t                       numPoloidalSlices,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r0,
               real_t                       r1,
               real_t                       r2 )
   : numToroidalSlices_( numToroidalSlices )
   , numPoloidalSlices_( numPoloidalSlices )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )
   , delta_( delta )
   , r0_( r0 )
   , r1_( r1 )
   , r2_( r2 )
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
      auto xold_0 = xold[0];
      auto xold_1 = xold[1];
      auto xold_2 = xold[2];
      auto tmp0   = pow( pow( xold_0, 2 ) + pow( xold_1, 2 ), -1.0 / 2.0 );
      auto tmp1   = radiusOriginToCenterOfTube_ * tmp0;
      auto tmp2   = 0.5 * toroidalAngleIncrement_;
      auto tmp3   = sin( tmp2 + toroidalAngleIncrement_ * toroidalPrism_ + toroidalStartAngle_ - atan2( xold_1, xold_0 ) +
                       1.5707963267948966 ) /
                  sin( tmp2 - 1.5707963267948966 );
      auto tmp4  = -tmp1 * xold_0 - tmp3 * xold_0;
      auto tmp5  = tmp0 * xold_0;
      auto tmp6  = -tmp1 * xold_1 - tmp3 * xold_1;
      auto tmp7  = tmp0 * xold_1;
      auto tmp8  = atan2( xold_2, tmp4 * tmp5 + tmp6 * tmp7 );
      auto tmp9  = ( ( tmp8 < 0 ) ? ( tmp8 + 6.2831853071795862 ) : ( tmp8 ) );
      auto tmp10 = sin( tmp9 );
      auto tmp11 = 0.5 * poloidalAngleIncrement_;
      auto tmp12 = sqrt( pow( tmp4, 2 ) + pow( tmp6, 2 ) + pow( xold_2, 2 ) ) *
                   sin( poloidalAngleIncrement_ * poloidalPrism_ + poloidalStartAngle_ + tmp11 - tmp9 + 1.5707963267948966 ) /
                   ( r0_ * tubeLayerRadiiBack_ * sin( tmp11 - 1.5707963267948966 ) );
      auto tmp13 = -r1_ * tmp12 * cos( tmp10 * asin( delta_ ) + tmp9 ) + 1;
      xnew[0]    = tmp13 * tmp5;
      xnew[1]    = tmp13 * tmp7;
      xnew[2]    = -r2_ * tmp10 * tmp12;
   }

   real_t evalDF( const Point3D& xold, Matrix3r& DF ) const final
   {
      auto xold_0 = xold[0];
      auto xold_1 = xold[1];
      auto xold_2 = xold[2];
      auto tmp0   = pow( xold_0, 2 );
      auto tmp1   = pow( xold_1, 2 );
      auto tmp2   = tmp0 + tmp1;
      auto tmp3   = pow( tmp2, -1.0 / 2.0 );
      auto tmp4   = radiusOriginToCenterOfTube_ * tmp3;
      auto tmp5   = 0.5 * toroidalAngleIncrement_;
      auto tmp6   = 1.0 / sin( tmp5 - 1.5707963267948966 );
      auto tmp7 =
          tmp5 + toroidalAngleIncrement_ * toroidalPrism_ + toroidalStartAngle_ - atan2( xold_1, xold_0 ) + 1.5707963267948966;
      auto tmp8  = tmp6 * sin( tmp7 );
      auto tmp9  = -tmp4 * xold_0 - tmp8 * xold_0;
      auto tmp10 = tmp3 * tmp9;
      auto tmp11 = -tmp4 * xold_1 - tmp8 * xold_1;
      auto tmp12 = tmp11 * tmp3;
      auto tmp13 = tmp10 * xold_0 + tmp12 * xold_1;
      auto tmp14 = atan2( xold_2, tmp13 );
      auto tmp15 = -poloidalStartAngle_ + tmp14;
      auto tmp16 = ( ( poloidalStartAngle_ - tmp14 > 0 ) ? ( tmp15 + 6.2831853071795862 ) : ( tmp15 ) );
      auto tmp17 = asin( delta_ );
      auto tmp18 = sin( tmp16 );
      auto tmp19 = tmp16 + tmp17 * tmp18;
      auto tmp20 = cos( tmp19 );
      auto tmp21 = pow( xold_2, 2 );
      auto tmp22 = sqrt( pow( tmp11, 2 ) + tmp21 + pow( tmp9, 2 ) );
      auto tmp23 = 1.0 / r0_;
      auto tmp24 = 1.0 / tubeLayerRadiiBack_;
      auto tmp25 = 0.5 * poloidalAngleIncrement_;
      auto tmp26 = 1.0 / sin( tmp25 - 1.5707963267948966 );
      auto tmp27 = poloidalAngleIncrement_ * poloidalPrism_ + poloidalStartAngle_ - tmp16 + tmp25 + 1.5707963267948966;
      auto tmp28 = tmp23 * tmp24 * tmp26 * sin( tmp27 );
      auto tmp29 = tmp22 * tmp28;
      auto tmp30 = r1_ * tmp29;
      auto tmp31 = -tmp20 * tmp30 + 1;
      auto tmp32 = tmp3 * tmp31;
      auto tmp33 = pow( tmp2, -3.0 / 2.0 );
      auto tmp34 = tmp0 * tmp33;
      auto tmp35 = r1_ * tmp20;
      auto tmp36 = xold_0 * xold_1;
      auto tmp37 = tmp33 * tmp36;
      auto tmp38 = radiusOriginToCenterOfTube_ * tmp37;
      auto tmp39 = 2 * tmp38;
      auto tmp40 = tmp6 * cos( tmp7 ) / tmp2;
      auto tmp41 = tmp1 * tmp40;
      auto tmp42 = ( 1.0 / 2.0 ) * tmp11;
      auto tmp43 = radiusOriginToCenterOfTube_ * tmp34;
      auto tmp44 = tmp36 * tmp40;
      auto tmp45 = 2 * tmp44;
      auto tmp46 = -2 * tmp4 - 2 * tmp8;
      auto tmp47 = ( 1.0 / 2.0 ) * tmp9;
      auto tmp48 = tmp28 / tmp22;
      auto tmp49 = tmp48 * ( tmp42 * ( tmp39 - 2 * tmp41 ) + tmp47 * ( 2 * tmp43 - tmp45 + tmp46 ) );
      auto tmp50 = tmp3 * xold_1;
      auto tmp51 = -tmp4 - tmp8;
      auto tmp52 = tmp3 * xold_0;
      auto tmp53 = 1.0 / ( pow( tmp13, 2 ) + tmp21 );
      auto tmp54 = tmp53 * xold_2;
      auto tmp55 =
          tmp54 * ( tmp10 - tmp11 * tmp37 - tmp34 * tmp9 + tmp50 * ( tmp38 - tmp41 ) + tmp52 * ( tmp43 - tmp44 + tmp51 ) );
      auto tmp56 = tmp22 * tmp23 * tmp24 * tmp26 * cos( tmp27 );
      auto tmp57 = tmp35 * tmp56;
      auto tmp58 = cos( tmp16 );
      auto tmp59 = tmp17 * tmp58;
      auto tmp60 = tmp30 * sin( tmp19 );
      auto tmp61 = -tmp35 * tmp49 - tmp55 * tmp57 + tmp60 * ( -tmp55 * tmp59 - tmp55 );
      auto tmp62 = -tmp31 * tmp37;
      auto tmp63 = tmp0 * tmp40;
      auto tmp64 = tmp1 * tmp33;
      auto tmp65 = radiusOriginToCenterOfTube_ * tmp64;
      auto tmp66 = tmp42 * ( tmp45 + tmp46 + 2 * tmp65 ) + tmp47 * ( tmp39 + 2 * tmp63 );
      auto tmp67 = tmp35 * tmp48;
      auto tmp68 =
          tmp54 * ( -tmp11 * tmp64 + tmp12 - tmp37 * tmp9 + tmp50 * ( tmp44 + tmp51 + tmp65 ) + tmp52 * ( tmp38 + tmp63 ) );
      auto tmp69 = -tmp57 * tmp68 + tmp60 * ( -tmp59 * tmp68 - tmp68 ) - tmp66 * tmp67;
      auto tmp70 = tmp13 * tmp53;
      auto tmp71 = tmp56 * tmp70;
      auto tmp72 = tmp35 * tmp71 + tmp60 * ( tmp59 * tmp70 + tmp70 ) - tmp67 * xold_2;
      auto tmp73 = r2_ * tmp18;
      auto tmp74 = r2_ * tmp29 * tmp58;
      auto tmp75 = tmp56 * tmp73;
      auto tmp76 = tmp48 * tmp73;
      DF( 0, 0 ) = -tmp31 * tmp34 + tmp32 + tmp52 * tmp61;
      DF( 0, 1 ) = tmp52 * tmp69 + tmp62;
      DF( 0, 2 ) = tmp52 * tmp72;
      DF( 1, 0 ) = tmp50 * tmp61 + tmp62;
      DF( 1, 1 ) = -tmp31 * tmp64 + tmp32 + tmp50 * tmp69;
      DF( 1, 2 ) = tmp50 * tmp72;
      DF( 2, 0 ) = -tmp49 * tmp73 + tmp55 * tmp74 - tmp55 * tmp75;
      DF( 2, 1 ) = -tmp66 * tmp76 + tmp68 * tmp74 - tmp68 * tmp75;
      DF( 2, 2 ) = -tmp70 * tmp74 + tmp71 * tmp73 - tmp76 * xold_2;
   }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const
   {
      WALBERLA_ABORT( "Serialization not implemented for TokamakMap" );
   }

   static void setMap( SetupPrimitiveStorage& setupStorage,
                       uint_t                 numToroidalSlices,
                       uint_t                 numPoloidalSlices,
                       real_t                 radiusOriginToCenterOfTube,
                       std::vector< real_t >  tubeLayerRadii,
                       real_t                 toroidalStartAngle,
                       real_t                 poloidalStartAngle,
                       real_t                 delta,
                       real_t                 r0,
                       real_t                 r1,
                       real_t                 r2 )
   {
      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(),
                                      std::make_shared< TokamakMap >( cell,
                                                                      setupStorage,
                                                                      numToroidalSlices,
                                                                      numPoloidalSlices,
                                                                      radiusOriginToCenterOfTube,
                                                                      tubeLayerRadii,
                                                                      toroidalStartAngle,
                                                                      poloidalStartAngle,
                                                                      delta,
                                                                      r0,
                                                                      r1,
                                                                      r2 ) );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(),
                                      std::make_shared< TokamakMap >( face,
                                                                      setupStorage,
                                                                      numToroidalSlices,
                                                                      numPoloidalSlices,
                                                                      radiusOriginToCenterOfTube,
                                                                      tubeLayerRadii,
                                                                      toroidalStartAngle,
                                                                      poloidalStartAngle,
                                                                      delta,
                                                                      r0,
                                                                      r1,
                                                                      r2 ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(),
                                      std::make_shared< TokamakMap >( edge,
                                                                      setupStorage,
                                                                      numToroidalSlices,
                                                                      numPoloidalSlices,
                                                                      radiusOriginToCenterOfTube,
                                                                      tubeLayerRadii,
                                                                      toroidalStartAngle,
                                                                      poloidalStartAngle,
                                                                      delta,
                                                                      r0,
                                                                      r1,
                                                                      r2 ) );
      }

      for ( auto it : setupStorage.getVertices() )
      {
         Vertex& vertex = *it.second;
         setupStorage.setGeometryMap( vertex.getID(),
                                      std::make_shared< TokamakMap >( vertex,
                                                                      setupStorage,
                                                                      numToroidalSlices,
                                                                      numPoloidalSlices,
                                                                      radiusOriginToCenterOfTube,
                                                                      tubeLayerRadii,
                                                                      toroidalStartAngle,
                                                                      poloidalStartAngle,
                                                                      delta,
                                                                      r0,
                                                                      r1,
                                                                      r2 ) );
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
   real_t angle( Point3D a, Point3D b ) const { return std::acos( a.dot( b ) / ( a.norm() * b.norm() ) ); }

   void identifyPrism( const Cell& cell )
   {
      tubeLayerRadiiBack_ = tubeLayerRadii_.back();

      toroidalAngleIncrement_ = 2 * pi / real_c( numToroidalSlices_ );
      poloidalAngleIncrement_ = 2 * pi / real_c( numPoloidalSlices_ );

      auto    coords = cell.getCoordinates();
      Point3D centroid( { 0, 0, 0 } );
      for ( uint_t i = 0; i < 4; i++ )
      {
         centroid += coords[i];
      }
      centroid *= 0.25;

      auto toroidalAngle = atan2( centroid[1], centroid[0] );
      toroidalAngle -= toroidalStartAngle_;
      if ( toroidalAngle < 0 )
      {
         toroidalAngle += 2 * pi;
      }
      toroidalPrism_ = uint_c( toroidalAngle / toroidalAngleIncrement_ );

      // to find out the poloidal prism we need to do some more work ...

      // first project the centroid to the torus (toroidal)
      {
         auto alpha = toroidalAngle - toroidalStartAngle_ - real_c( toroidalPrism_ ) * toroidalAngleIncrement_;
         auto beta  = 0.5 * ( pi - toroidalAngleIncrement_ );
         auto gamma = pi - alpha - beta;
         auto toroidalRadiusNew =
             ( std::sin( gamma ) * ( std::sqrt( centroid[0] * centroid[0] + centroid[1] * centroid[1] ) / std::sin( beta ) ) );

         centroid[0] = toroidalRadiusNew * std::cos( toroidalAngle );
         centroid[1] = toroidalRadiusNew * std::sin( toroidalAngle );
         centroid[2] = centroid[2];
      }

      // then we rotate the mapped centroid around the z-axis and translate it to the origin
      // this way we can find the angle and therefore the prism ID via polar coordinates in the x-z-plane

      auto C                     = torusCoordinates( radiusOriginToCenterOfTube_, 0, toroidalAngle, 0 );
      auto centroidTrafoToOrigin = centroid - C;
      centroidTrafoToOrigin      = Point3D(
          { std::cos( -toroidalAngle ) * centroidTrafoToOrigin[0] - std::sin( -toroidalAngle ) * centroidTrafoToOrigin[1],
            std::sin( -toroidalAngle ) * centroidTrafoToOrigin[0] + std::cos( -toroidalAngle ) * centroidTrafoToOrigin[1],
            centroidTrafoToOrigin[2] } );

      auto poloidalAngle = std::atan2( centroidTrafoToOrigin[2], centroidTrafoToOrigin[0] );
      poloidalAngle -= poloidalStartAngle_;
      if ( poloidalAngle < 0 )
      {
         poloidalAngle += 2 * pi;
      }
      poloidalPrism_ = uint_c( ( poloidalAngle ) / poloidalAngleIncrement_ );
   }

   uint_t                numToroidalSlices_;
   uint_t                numPoloidalSlices_;
   real_t                radiusOriginToCenterOfTube_;
   std::vector< real_t > tubeLayerRadii_;
   real_t                toroidalStartAngle_;
   real_t                poloidalStartAngle_;

   real_t toroidalAngleIncrement_;
   real_t poloidalAngleIncrement_;

   uint_t toroidalPrism_;
   uint_t poloidalPrism_;

   Point3D sliceCenterFront_;
   Point3D sliceCenterBack_;

   real_t delta_;
   real_t r0_;
   real_t r1_;
   real_t r2_;

   real_t tubeLayerRadiiBack_;
};

} // end of namespace hyteg
