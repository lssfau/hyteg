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
#include <utility>

#include "core/DataTypes.h"
#include "core/math/Constants.h"

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/geometry/Polygons.hpp"
#include "hyteg/geometry/Torus.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::int_c;
using walberla::real_c;
using walberla::math::pi;

namespace hyteg {

/// The class implements a mapping of the torus mesh to a tokamak domain.
class TokamakMap : public GeometryMap
{
 public:
   TokamakMap( const Cell& cell,
               const SetupPrimitiveStorage&,
               uint_t                toroidalResolution,
               uint_t                poloidalResolution,
               real_t                radiusOriginToCenterOfTube,
               std::vector< real_t > tubeLayerRadii,
               real_t                toroidalStartAngle,
               real_t                poloidalStartAngle,
               real_t                delta,
               real_t                r1,
               real_t                r2 )
   : toroidalResolution_( toroidalResolution )
   , poloidalResolution_( poloidalResolution )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( std::move( tubeLayerRadii ) )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )
   , delta_( delta )
   , r1_( r1 )
   , r2_( r2 )

   {
      identifyPrism( cell );
   }

   TokamakMap( const Face&                  face,
               const SetupPrimitiveStorage& setupStorage,
               uint_t                       toroidalResolution,
               uint_t                       poloidalResolution,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r1,
               real_t                       r2 )
   : toroidalResolution_( toroidalResolution )
   , poloidalResolution_( poloidalResolution )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( std::move( tubeLayerRadii ) )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )
   , delta_( delta )
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
               uint_t                       toroidalResolution,
               uint_t                       poloidalResolution,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r1,
               real_t                       r2 )
   : toroidalResolution_( toroidalResolution )
   , poloidalResolution_( poloidalResolution )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( std::move( tubeLayerRadii ) )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )
   , delta_( delta )
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
               uint_t                       toroidalResolution,
               uint_t                       poloidalResolution,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r1,
               real_t                       r2 )
   : toroidalResolution_( toroidalResolution )
   , poloidalResolution_( poloidalResolution )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( std::move( tubeLayerRadii ) )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )
   , delta_( delta )
   , r1_( r1 )
   , r2_( r2 )
   {
      std::vector< PrimitiveID > neighborCells;
      vertex.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TokamakMap( walberla::mpi::RecvBuffer& recvBuffer )
   {
      recvBuffer >> toroidalResolution_;
      recvBuffer >> poloidalResolution_;
      recvBuffer >> radiusOriginToCenterOfTube_;
      recvBuffer >> tubeLayerRadii_;
      recvBuffer >> toroidalStartAngle_;
      recvBuffer >> poloidalStartAngle_;

      recvBuffer >> toroidalAngleIncrement_;
      recvBuffer >> poloidalAngleIncrement_;

      recvBuffer >> toroidalPrism_;
      recvBuffer >> poloidalPrism_;

      recvBuffer >> sliceCenterFront_;
      recvBuffer >> sliceCenterBack_;

      recvBuffer >> delta_;
      recvBuffer >> r1_;
      recvBuffer >> r2_;

      recvBuffer >> tubeLayerRadiiBack_;
   }

   void evalF( const Point3D& xold, Point3D& xnew ) const final
   {
      // generated with data/codegen/geometry/TokamakMap.py

      auto xold_0 = xold[0];
      auto xold_1 = xold[1];
      auto xold_2 = xold[2];
      auto tmp0   = pow( xold_2, 2 );
      auto tmp1   = pow( pow( xold_0, 2 ) + pow( xold_1, 2 ), -1.0 / 2.0 );
      auto tmp2   = tmp1 * xold_0;
      auto tmp3   = 0.5 * toroidalAngleIncrement_;
      auto tmp4 = sin( tmp3 + toroidalAngleIncrement_ * real_c( toroidalPrism_ ) + toroidalStartAngle_ - atan2( xold_1, xold_0 ) +
                       1.5707963267948966 ) /
                  sin( tmp3 - 1.5707963267948966 );
      auto tmp5  = -radiusOriginToCenterOfTube_ * tmp2 - tmp4 * xold_0;
      auto tmp6  = tmp1 * xold_1;
      auto tmp7  = -radiusOriginToCenterOfTube_ * tmp6 - tmp4 * xold_1;
      auto tmp8  = tmp2 * tmp5 + tmp6 * tmp7;
      auto tmp9  = sqrt( tmp0 + pow( tmp8, 2 ) ) < 1.0e-14;
      auto tmp10 = atan2( xold_2, tmp8 );
      auto tmp11 = -poloidalStartAngle_ + tmp10;
      auto tmp12 = ( ( poloidalStartAngle_ - tmp10 > 0 ) ? ( tmp11 + 6.2831853071795862 ) : ( tmp11 ) );
      auto tmp13 = poloidalStartAngle_ + ( ( tmp9 ) ? ( 0 ) : ( tmp12 ) );
      auto tmp14 = sin( tmp13 );
      auto tmp15 = 0.5 * poloidalAngleIncrement_;
      auto tmp16 = poloidalAngleIncrement_ * real_c( poloidalPrism_ );
      auto tmp17 =
          ( ( tmp9 ) ?
                ( 0 ) :
                ( -sqrt( tmp0 + pow( tmp5, 2 ) + pow( tmp7, 2 ) ) *
                  sin( -tmp12 + tmp15 + ( ( tmp16 < 0 ) ? ( tmp16 + 6.2831853071795862 ) : ( tmp16 ) ) + 1.5707963267948966 ) /
                  sin( tmp15 - 1.5707963267948966 ) ) ) /
          tubeLayerRadiiBack_;
      auto tmp18 = r1_ * tmp17 * cos( tmp13 + tmp14 * asin( delta_ ) ) + radiusOriginToCenterOfTube_;
      xnew[0]    = real_c( tmp18 * tmp2 );
      xnew[1]    = real_c( tmp18 * tmp6 );
      xnew[2]    = real_c( r2_ * tmp14 * tmp17 );
   }

   real_t evalDF( const Point3D& xold, Matrix3r& DF ) const final
   {
      // generated with data/codegen/geometry/TokamakMap.py

      auto xold_0 = xold[0];
      auto xold_1 = xold[1];
      auto xold_2 = xold[2];
      auto tmp0   = pow( xold_0, 2 );
      auto tmp1   = pow( xold_1, 2 );
      auto tmp2   = tmp0 + tmp1;
      auto tmp3   = pow( tmp2, -1.0 / 2.0 );
      auto tmp4   = 1.0 / tubeLayerRadiiBack_;
      auto tmp5   = pow( xold_2, 2 );
      auto tmp6   = radiusOriginToCenterOfTube_ * tmp3;
      auto tmp7   = 0.5 * toroidalAngleIncrement_;
      auto tmp8   = 1.0 / sin( tmp7 - 1.5707963267948966 );
      auto tmp9   = tmp7 + toroidalAngleIncrement_ * real_c( toroidalPrism_ ) + toroidalStartAngle_ - atan2( xold_1, xold_0 ) +
                  1.5707963267948966;
      auto tmp10 = tmp8 * sin( tmp9 );
      auto tmp11 = -tmp10 * xold_0 - tmp6 * xold_0;
      auto tmp12 = tmp11 * tmp3;
      auto tmp13 = -tmp10 * xold_1 - tmp6 * xold_1;
      auto tmp14 = tmp13 * tmp3;
      auto tmp15 = tmp12 * xold_0 + tmp14 * xold_1;
      auto tmp16 = pow( tmp15, 2 ) + tmp5;
      auto tmp17 = sqrt( tmp16 ) < 1.0e-14;
      auto tmp18 = sqrt( pow( tmp11, 2 ) + pow( tmp13, 2 ) + tmp5 );
      auto tmp19 = 0.5 * poloidalAngleIncrement_;
      auto tmp20 = 1.0 / sin( tmp19 - 1.5707963267948966 );
      auto tmp21 = poloidalAngleIncrement_ * real_c( poloidalPrism_ );
      auto tmp22 = atan2( xold_2, tmp15 );
      auto tmp23 = -poloidalStartAngle_ + tmp22;
      auto tmp24 = ( ( poloidalStartAngle_ - tmp22 > 0 ) ? ( tmp23 + 6.2831853071795862 ) : ( tmp23 ) );
      auto tmp25 = tmp19 - tmp24 + ( ( tmp21 < 0 ) ? ( tmp21 + 6.2831853071795862 ) : ( tmp21 ) ) + 1.5707963267948966;
      auto tmp26 = tmp20 * sin( tmp25 );
      auto tmp27 = tmp4 * ( ( tmp17 ) ? ( 0 ) : ( -tmp18 * tmp26 ) );
      auto tmp28 = asin( delta_ );
      auto tmp29 = poloidalStartAngle_ + ( ( tmp17 ) ? ( 0 ) : ( tmp24 ) );
      auto tmp30 = sin( tmp29 );
      auto tmp31 = tmp28 * tmp30 + tmp29;
      auto tmp32 = r1_ * cos( tmp31 );
      auto tmp33 = radiusOriginToCenterOfTube_ + tmp27 * tmp32;
      auto tmp34 = tmp3 * tmp33;
      auto tmp35 = pow( tmp2, -3.0 / 2.0 );
      auto tmp36 = tmp0 * tmp35;
      auto tmp37 = xold_0 * xold_1;
      auto tmp38 = tmp35 * tmp37;
      auto tmp39 = radiusOriginToCenterOfTube_ * tmp38;
      auto tmp40 = 2 * tmp39;
      auto tmp41 = tmp8 * cos( tmp9 ) / tmp2;
      auto tmp42 = tmp1 * tmp41;
      auto tmp43 = ( 1.0 / 2.0 ) * tmp13;
      auto tmp44 = radiusOriginToCenterOfTube_ * tmp36;
      auto tmp45 = tmp37 * tmp41;
      auto tmp46 = 2 * tmp45;
      auto tmp47 = -2 * tmp10 - 2 * tmp6;
      auto tmp48 = ( 1.0 / 2.0 ) * tmp11;
      auto tmp49 = tmp26 / tmp18;
      auto tmp50 = tmp3 * xold_1;
      auto tmp51 = -tmp10 - tmp6;
      auto tmp52 = tmp3 * xold_0;
      auto tmp53 = 1.0 / tmp16;
      auto tmp54 = tmp53 * xold_2;
      auto tmp55 =
          tmp54 * ( -tmp11 * tmp36 + tmp12 - tmp13 * tmp38 + tmp50 * ( tmp39 - tmp42 ) + tmp52 * ( tmp44 - tmp45 + tmp51 ) );
      auto tmp56 = tmp18 * tmp20 * cos( tmp25 );
      auto tmp57 =
          ( ( tmp17 ) ? ( 0 ) :
                        ( -tmp49 * ( tmp43 * ( tmp40 - 2 * tmp42 ) + tmp48 * ( 2 * tmp44 - tmp46 + tmp47 ) ) - tmp55 * tmp56 ) );
      auto tmp58 = tmp32 * tmp4;
      auto tmp59 = ( ( tmp17 ) ? ( 0 ) : ( -tmp55 ) );
      auto tmp60 = cos( tmp29 );
      auto tmp61 = tmp28 * tmp60;
      auto tmp62 = r1_ * tmp27 * sin( tmp31 );
      auto tmp63 = tmp57 * tmp58 - tmp62 * ( tmp59 * tmp61 + tmp59 );
      auto tmp64 = -tmp33 * tmp38;
      auto tmp65 = tmp0 * tmp41;
      auto tmp66 = tmp1 * tmp35;
      auto tmp67 = radiusOriginToCenterOfTube_ * tmp66;
      auto tmp68 =
          tmp54 * ( -tmp11 * tmp38 - tmp13 * tmp66 + tmp14 + tmp50 * ( tmp45 + tmp51 + tmp67 ) + tmp52 * ( tmp39 + tmp65 ) );
      auto tmp69 =
          ( ( tmp17 ) ? ( 0 ) :
                        ( -tmp49 * ( tmp43 * ( tmp46 + tmp47 + 2 * tmp67 ) + tmp48 * ( tmp40 + 2 * tmp65 ) ) - tmp56 * tmp68 ) );
      auto tmp70 = ( ( tmp17 ) ? ( 0 ) : ( -tmp68 ) );
      auto tmp71 = tmp58 * tmp69 - tmp62 * ( tmp61 * tmp70 + tmp70 );
      auto tmp72 = tmp15 * tmp53;
      auto tmp73 = ( ( tmp17 ) ? ( 0 ) : ( -tmp49 * xold_2 + tmp56 * tmp72 ) );
      auto tmp74 = ( ( tmp17 ) ? ( 0 ) : ( tmp72 ) );
      auto tmp75 = tmp58 * tmp73 - tmp62 * ( tmp61 * tmp74 + tmp74 );
      auto tmp76 = r2_ * tmp27 * tmp60;
      auto tmp77 = r2_ * tmp30 * tmp4;
      DF( 0, 0 ) = real_c( -tmp33 * tmp36 + tmp34 + tmp52 * tmp63 );
      DF( 0, 1 ) = real_c( tmp52 * tmp71 + tmp64 );
      DF( 0, 2 ) = real_c( tmp52 * tmp75 );
      DF( 1, 0 ) = real_c( tmp50 * tmp63 + tmp64 );
      DF( 1, 1 ) = real_c( -tmp33 * tmp66 + tmp34 + tmp50 * tmp71 );
      DF( 1, 2 ) = real_c( tmp50 * tmp75 );
      DF( 2, 0 ) = real_c( tmp57 * tmp77 + tmp59 * tmp76 );
      DF( 2, 1 ) = real_c( tmp69 * tmp77 + tmp70 * tmp76 );
      DF( 2, 2 ) = real_c( tmp73 * tmp77 + tmp74 * tmp76 );

      WALBERLA_ASSERT( !std::isnan( DF( 0, 0 ) ), "The tokamap map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 0, 1 ) ), "The tokamap map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 0, 2 ) ), "The tokamap map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 1, 0 ) ), "The tokamap map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 1, 1 ) ), "The tokamap map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 1, 2 ) ), "The tokamap map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 2, 0 ) ), "The tokamap map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 2, 1 ) ), "The tokamap map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 2, 2 ) ), "The tokamap map produces NaNs :( again." )

      return DF.determinant();
   }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const final
   {
      sendBuffer << Type::TOKAMAK;
      sendBuffer << toroidalResolution_;
      sendBuffer << poloidalResolution_;
      sendBuffer << radiusOriginToCenterOfTube_;
      sendBuffer << tubeLayerRadii_;
      sendBuffer << toroidalStartAngle_;
      sendBuffer << poloidalStartAngle_;

      sendBuffer << toroidalAngleIncrement_;
      sendBuffer << poloidalAngleIncrement_;

      sendBuffer << toroidalPrism_;
      sendBuffer << poloidalPrism_;

      sendBuffer << sliceCenterFront_;
      sendBuffer << sliceCenterBack_;

      sendBuffer << delta_;
      sendBuffer << r1_;
      sendBuffer << r2_;

      sendBuffer << tubeLayerRadiiBack_;
   }

   /// \brief Applies the tokamak map to the SetupPrimitiveStorage.
   ///
   /// As unstructured base mesh, please use the mesh generator MeshInfo::meshTorus() with the same parameters for the torus.
   ///
   /// The tokamak geometry description is taken from
   ///
   ///     Oueslati et al. (2019) Numerical derivation of steady flows in visco-resistive magnetohydrodynamics
   ///                            for JET and ITER-like geometries with no symmetry breaking.
   ///
   /// Can be abused to map to a standard torus by setting delta == 0 and r1 == r2 == tubeLayerRadii.back()
   ///
   ///     >>> THERE IS A DEDICATED TORUS MAP THAT AVOIDS UNNECESSARY FLOPS, THOUGH! <<<<
   ///
   /// \param setupStorage the SetupPrimitiveStorage instance
   /// \param toroidalResolution number of prisms in toroidal direction (along the ring) in a complete (360 degree) ring
   /// \param poloidalResolution number of vertices on the boundary of a slice through the tube
   /// \param radiusOriginToCenterOfTube distance from origin to the center of the tube
   /// \param tubeLayerRadii list of radii of layers of the sliced tube - the last element defines the actual radius of the tube
   /// \param toroidalStartAngle angle (in radians) by which the domain shall be rotated about the z-axis
   /// \param poloidalStartAngle angle (in radians) by which the domain shall be rotated about the ring through the center of the tube
   /// \param delta triangularity parameter of the tokamak
   /// \param r1 semi-minor axis radius
   /// \param r2 semi-major axis radius
   ///
   static void setMap( SetupPrimitiveStorage&       setupStorage,
                       uint_t                       toroidalResolution,
                       uint_t                       poloidalResolution,
                       real_t                       radiusOriginToCenterOfTube,
                       const std::vector< real_t >& tubeLayerRadii,
                       real_t                       toroidalStartAngle,
                       real_t                       poloidalStartAngle,
                       real_t                       delta,
                       real_t                       r1,
                       real_t                       r2 )
   {
      for ( const auto& it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(),
                                      std::make_shared< TokamakMap >( cell,
                                                                      setupStorage,
                                                                      toroidalResolution,
                                                                      poloidalResolution,
                                                                      radiusOriginToCenterOfTube,
                                                                      tubeLayerRadii,
                                                                      toroidalStartAngle,
                                                                      poloidalStartAngle,
                                                                      delta,
                                                                      r1,
                                                                      r2 ) );
      }

      for ( const auto& it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(),
                                      std::make_shared< TokamakMap >( face,
                                                                      setupStorage,
                                                                      toroidalResolution,
                                                                      poloidalResolution,
                                                                      radiusOriginToCenterOfTube,
                                                                      tubeLayerRadii,
                                                                      toroidalStartAngle,
                                                                      poloidalStartAngle,
                                                                      delta,
                                                                      r1,
                                                                      r2 ) );
      }

      for ( const auto& it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(),
                                      std::make_shared< TokamakMap >( edge,
                                                                      setupStorage,
                                                                      toroidalResolution,
                                                                      poloidalResolution,
                                                                      radiusOriginToCenterOfTube,
                                                                      tubeLayerRadii,
                                                                      toroidalStartAngle,
                                                                      poloidalStartAngle,
                                                                      delta,
                                                                      r1,
                                                                      r2 ) );
      }

      for ( const auto& it : setupStorage.getVertices() )
      {
         Vertex& vertex = *it.second;
         setupStorage.setGeometryMap( vertex.getID(),
                                      std::make_shared< TokamakMap >( vertex,
                                                                      setupStorage,
                                                                      toroidalResolution,
                                                                      poloidalResolution,
                                                                      radiusOriginToCenterOfTube,
                                                                      tubeLayerRadii,
                                                                      toroidalStartAngle,
                                                                      poloidalStartAngle,
                                                                      delta,
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
   static real_t angle( const Point3D& a, const Point3D& b ) { return std::acos( a.dot( b ) / ( a.norm() * b.norm() ) ); }

   void identifyPrism( const Cell& cell )
   {
      tubeLayerRadiiBack_ = tubeLayerRadii_.back();

      toroidalAngleIncrement_ = 2 * pi / real_c( toroidalResolution_ );
      poloidalAngleIncrement_ = 2 * pi / real_c( poloidalResolution_ );

      auto    coords = cell.getCoordinates();
      Point3D centroid( 0, 0, 0 );
      for ( uint_t i = 0; i < 4; i++ )
      {
         centroid += coords[i];
      }
      centroid *= 0.25;

      real_t toroidalAngle = real_c( atan2( centroid[1], centroid[0] ) );
      toroidalAngle -= toroidalStartAngle_;
      if ( toroidalAngle < 0 )
      {
         toroidalAngle += 2 * pi;
      }
      toroidalPrism_ = uint_c( toroidalAngle / toroidalAngleIncrement_ );

      // to find out the poloidal prism we need to do some more work ...

      // first project the centroid to the torus (toroidal)
      {
         real_t alpha             = toroidalAngle - toroidalStartAngle_ - real_c( toroidalPrism_ ) * toroidalAngleIncrement_;
         real_t beta              = real_c( 0.5 ) * ( pi - toroidalAngleIncrement_ );
         real_t gamma             = pi - alpha - beta;
         real_t toroidalRadiusNew = real_c(
             std::sin( gamma ) * ( std::sqrt( centroid[0] * centroid[0] + centroid[1] * centroid[1] ) / std::sin( beta ) ) );

         centroid[0] = toroidalRadiusNew * real_c( std::cos( toroidalAngle ) );
         centroid[1] = toroidalRadiusNew * real_c( std::sin( toroidalAngle ) );
         centroid[2] = centroid[2];
      }

      // then we rotate the mapped centroid around the z-axis and translate it to the origin
      // this way we can find the angle and therefore the prism ID via polar coordinates in the x-z-plane

      auto    C                     = torusCoordinates( radiusOriginToCenterOfTube_, 0, toroidalAngle, 0 );
      Point3D centroidTrafoToOrigin = centroid - C;
      centroidTrafoToOrigin         = Point3D( { real_c( std::cos( -toroidalAngle ) * centroidTrafoToOrigin[0] -
                                                 std::sin( -toroidalAngle ) * centroidTrafoToOrigin[1] ),
                                                 real_c( std::sin( -toroidalAngle ) * centroidTrafoToOrigin[0] +
                                                 std::cos( -toroidalAngle ) * centroidTrafoToOrigin[1] ),
                                                 centroidTrafoToOrigin[2] } );

      auto poloidalAngle = std::atan2( centroidTrafoToOrigin[2], centroidTrafoToOrigin[0] );
      poloidalAngle -= poloidalStartAngle_;
      if ( poloidalAngle < 0 )
      {
         poloidalAngle += 2 * pi;
      }
      poloidalPrism_ = uint_c( ( poloidalAngle ) / poloidalAngleIncrement_ );
   }

   uint_t                toroidalResolution_;
   uint_t                poloidalResolution_;
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
   real_t r1_;
   real_t r2_;

   real_t tubeLayerRadiiBack_;
};

} // end of namespace hyteg
