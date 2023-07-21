/*
* Copyright (c) 2017-2023 Nils Kohl.
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

/// The class implements a mapping of the torus mesh to a (curved) torus domain.
class TorusMap : public GeometryMap
{
 public:
   TorusMap( const Cell&                  cell,
             const SetupPrimitiveStorage& setupStorage,
             uint_t                       toroidalResolution,
             uint_t                       poloidalResolution,
             real_t                       radiusOriginToCenterOfTube,
             std::vector< real_t >        tubeLayerRadii,
             real_t                       toroidalStartAngle,
             real_t                       poloidalStartAngle )
   : toroidalResolution( toroidalResolution )
   , poloidalResolution( poloidalResolution )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )

   {
      identifyPrism( cell );
   }

   TorusMap( const Face&                  face,
             const SetupPrimitiveStorage& setupStorage,
             uint_t                       toroidalResolution,
             uint_t                       poloidalResolution,
             real_t                       radiusOriginToCenterOfTube,
             std::vector< real_t >        tubeLayerRadii,
             real_t                       toroidalStartAngle,
             real_t                       poloidalStartAngle )
   : toroidalResolution( toroidalResolution )
   , poloidalResolution( poloidalResolution )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )

   {
      std::vector< PrimitiveID > neighborCells;
      face.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TorusMap( const Edge&                  edge,
             const SetupPrimitiveStorage& setupStorage,
             uint_t                       toroidalResolution,
             uint_t                       poloidalResolution,
             real_t                       radiusOriginToCenterOfTube,
             std::vector< real_t >        tubeLayerRadii,
             real_t                       toroidalStartAngle,
             real_t                       poloidalStartAngle )
   : toroidalResolution( toroidalResolution )
   , poloidalResolution( poloidalResolution )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )

   {
      std::vector< PrimitiveID > neighborCells;
      edge.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TorusMap( const Vertex&                vertex,
             const SetupPrimitiveStorage& setupStorage,
             uint_t                       toroidalResolution,
             uint_t                       poloidalResolution,
             real_t                       radiusOriginToCenterOfTube,
             std::vector< real_t >        tubeLayerRadii,
             real_t                       toroidalStartAngle,
             real_t                       poloidalStartAngle )
   : toroidalResolution( toroidalResolution )
   , poloidalResolution( poloidalResolution )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
   , toroidalStartAngle_( toroidalStartAngle )
   , poloidalStartAngle_( poloidalStartAngle )

   {
      std::vector< PrimitiveID > neighborCells;
      vertex.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *setupStorage.getCell( neighborCells[0] );
      identifyPrism( cell );
   }

   TorusMap( walberla::mpi::RecvBuffer& recvBuffer ) { WALBERLA_ABORT( "Deserialization not implemented for TorusMap" ); }

   void evalF( const Point3D& xold, Point3D& xnew ) const override final
   {
      // generated with data/codegen/geometry/TokamakMap.py
      // you need to specify that you want to skip the tokamak part inside the python file!

      auto xold_0 = xold[0];
      auto xold_1 = xold[1];
      auto xold_2 = xold[2];
      auto tmp0   = 0.5 * toroidalAngleIncrement_;
      auto tmp1   = sin( tmp0 + toroidalAngleIncrement_ * toroidalPrism_ + toroidalStartAngle_ - atan2( xold_1, xold_0 ) +
                       1.5707963267948966 ) /
                  sin( tmp0 - 1.5707963267948966 );
      auto tmp2  = -tmp1 * xold_0;
      auto tmp3  = pow( xold_2, 2 );
      auto tmp4  = pow( pow( xold_0, 2 ) + pow( xold_1, 2 ), -1.0 / 2.0 );
      auto tmp5  = tmp4 * xold_0;
      auto tmp6  = -radiusOriginToCenterOfTube_ * tmp5 + tmp2;
      auto tmp7  = tmp4 * xold_1;
      auto tmp8  = -tmp1 * xold_1;
      auto tmp9  = -radiusOriginToCenterOfTube_ * tmp7 + tmp8;
      auto tmp10 = tmp5 * tmp6 + tmp7 * tmp9;
      auto tmp11 = sqrt( pow( tmp10, 2 ) + tmp3 ) < 1.0e-14;
      auto tmp12 = atan2( xold_2, tmp10 );
      auto tmp13 = -poloidalStartAngle_ + tmp12;
      auto tmp14 = ( ( poloidalStartAngle_ - tmp12 > 0 ) ? ( tmp13 + 6.2831853071795862 ) : ( tmp13 ) );
      auto tmp15 = poloidalStartAngle_ + tmp14;
      auto tmp16 = 0.5 * poloidalAngleIncrement_;
      auto tmp17 = poloidalAngleIncrement_ * poloidalPrism_;
      auto tmp18 = sqrt( tmp3 + pow( tmp6, 2 ) + pow( tmp9, 2 ) ) *
                   sin( -tmp14 + tmp16 + ( ( tmp17 < 0 ) ? ( tmp17 + 6.2831853071795862 ) : ( tmp17 ) ) + 1.5707963267948966 ) /
                   sin( tmp16 - 1.5707963267948966 );
      auto tmp19 = radiusOriginToCenterOfTube_ - tmp18 * cos( tmp15 );
      xnew[0]    = ( ( tmp11 ) ? ( tmp2 ) : ( tmp19 * tmp5 ) );
      xnew[1]    = ( ( tmp11 ) ? ( tmp8 ) : ( tmp19 * tmp7 ) );
      xnew[2]    = ( ( tmp11 ) ? ( xold_2 ) : ( -tmp18 * sin( tmp15 ) ) );
   }

   real_t evalDF( const Point3D& xold, Matrix3r& DF ) const override final
   {
      // generated with data/codegen/geometry/TokamakMap.py
      // you need to specify that you want to skip the tokamak part inside the python file!

      auto xold_0 = xold[0];
      auto xold_1 = xold[1];
      auto xold_2 = xold[2];
      auto tmp0   = 0.5 * toroidalAngleIncrement_;
      auto tmp1   = 1.0 / sin( tmp0 - 1.5707963267948966 );
      auto tmp2 =
          tmp0 + toroidalAngleIncrement_ * toroidalPrism_ + toroidalStartAngle_ - atan2( xold_1, xold_0 ) + 1.5707963267948966;
      auto tmp3  = tmp1 * sin( tmp2 );
      auto tmp4  = -tmp3;
      auto tmp5  = xold_0 * xold_1;
      auto tmp6  = pow( xold_0, 2 );
      auto tmp7  = pow( xold_1, 2 );
      auto tmp8  = tmp6 + tmp7;
      auto tmp9  = tmp1 * cos( tmp2 ) / tmp8;
      auto tmp10 = tmp5 * tmp9;
      auto tmp11 = -tmp10 + tmp4;
      auto tmp12 = pow( xold_2, 2 );
      auto tmp13 = pow( tmp8, -1.0 / 2.0 );
      auto tmp14 = radiusOriginToCenterOfTube_ * tmp13;
      auto tmp15 = -tmp14 * xold_0 - tmp3 * xold_0;
      auto tmp16 = tmp13 * tmp15;
      auto tmp17 = -tmp14 * xold_1 - tmp3 * xold_1;
      auto tmp18 = tmp13 * tmp17;
      auto tmp19 = tmp16 * xold_0 + tmp18 * xold_1;
      auto tmp20 = tmp12 + pow( tmp19, 2 );
      auto tmp21 = sqrt( tmp20 ) < 1.0e-14;
      auto tmp22 = sqrt( tmp12 + pow( tmp15, 2 ) + pow( tmp17, 2 ) );
      auto tmp23 = atan2( xold_2, tmp19 );
      auto tmp24 = -poloidalStartAngle_ + tmp23;
      auto tmp25 = ( ( poloidalStartAngle_ - tmp23 > 0 ) ? ( tmp24 + 6.2831853071795862 ) : ( tmp24 ) );
      auto tmp26 = poloidalStartAngle_ + tmp25;
      auto tmp27 = cos( tmp26 );
      auto tmp28 = 0.5 * poloidalAngleIncrement_;
      auto tmp29 = 1.0 / sin( tmp28 - 1.5707963267948966 );
      auto tmp30 = poloidalAngleIncrement_ * poloidalPrism_;
      auto tmp31 = -tmp25 + tmp28 + ( ( tmp30 < 0 ) ? ( tmp30 + 6.2831853071795862 ) : ( tmp30 ) ) + 1.5707963267948966;
      auto tmp32 = tmp29 * sin( tmp31 );
      auto tmp33 = tmp27 * tmp32;
      auto tmp34 = tmp22 * tmp33;
      auto tmp35 = radiusOriginToCenterOfTube_ - tmp34;
      auto tmp36 = tmp13 * tmp35;
      auto tmp37 = pow( tmp8, -3.0 / 2.0 );
      auto tmp38 = tmp37 * tmp6;
      auto tmp39 = tmp37 * tmp5;
      auto tmp40 = radiusOriginToCenterOfTube_ * tmp39;
      auto tmp41 = 2 * tmp40;
      auto tmp42 = tmp7 * tmp9;
      auto tmp43 = ( 1.0 / 2.0 ) * tmp17;
      auto tmp44 = radiusOriginToCenterOfTube_ * tmp38;
      auto tmp45 = 2 * tmp10;
      auto tmp46 = -2 * tmp14 - 2 * tmp3;
      auto tmp47 = ( 1.0 / 2.0 ) * tmp15;
      auto tmp48 = tmp43 * ( tmp41 - 2 * tmp42 ) + tmp47 * ( 2 * tmp44 - tmp45 + tmp46 );
      auto tmp49 = 1.0 / tmp22;
      auto tmp50 = tmp33 * tmp49;
      auto tmp51 = sin( tmp26 );
      auto tmp52 = tmp32 * tmp51;
      auto tmp53 = -tmp42;
      auto tmp54 = tmp13 * xold_1;
      auto tmp55 = -tmp14;
      auto tmp56 = tmp13 * xold_0;
      auto tmp57 = -tmp15 * tmp38 + tmp16 - tmp17 * tmp39 + tmp54 * ( tmp40 + tmp53 ) + tmp56 * ( tmp11 + tmp44 + tmp55 );
      auto tmp58 = 1.0 / tmp20;
      auto tmp59 = tmp58 * xold_2;
      auto tmp60 = tmp22 * tmp59;
      auto tmp61 = tmp57 * tmp60;
      auto tmp62 = tmp29 * cos( tmp31 );
      auto tmp63 = tmp27 * tmp62;
      auto tmp64 = -tmp48 * tmp50 - tmp52 * tmp61 - tmp61 * tmp63;
      auto tmp65 = tmp6 * tmp9;
      auto tmp66 = -tmp35 * tmp39;
      auto tmp67 = tmp37 * tmp7;
      auto tmp68 = radiusOriginToCenterOfTube_ * tmp67;
      auto tmp69 = tmp43 * ( tmp45 + tmp46 + 2 * tmp68 ) + tmp47 * ( tmp41 + 2 * tmp65 );
      auto tmp70 = tmp10 + tmp4;
      auto tmp71 = -tmp15 * tmp39 - tmp17 * tmp67 + tmp18 + tmp54 * ( tmp55 + tmp68 + tmp70 ) + tmp56 * ( tmp40 + tmp65 );
      auto tmp72 = tmp60 * tmp71;
      auto tmp73 = -tmp50 * tmp69 - tmp52 * tmp72 - tmp63 * tmp72;
      auto tmp74 = tmp19 * tmp58;
      auto tmp75 = tmp22 * tmp74;
      auto tmp76 = tmp62 * tmp75;
      auto tmp77 = tmp27 * tmp76 - tmp50 * xold_2 + tmp52 * tmp75;
      auto tmp78 = tmp49 * tmp52;
      auto tmp79 = tmp34 * tmp59;
      auto tmp80 = tmp51 * tmp62;
      DF( 0, 0 ) = ( ( tmp21 ) ? ( tmp11 ) : ( -tmp35 * tmp38 + tmp36 + tmp56 * tmp64 ) );
      DF( 0, 1 ) = ( ( tmp21 ) ? ( tmp65 ) : ( tmp56 * tmp73 + tmp66 ) );
      DF( 0, 2 ) = ( ( tmp21 ) ? ( 0 ) : ( tmp56 * tmp77 ) );
      DF( 1, 0 ) = ( ( tmp21 ) ? ( tmp53 ) : ( tmp54 * tmp64 + tmp66 ) );
      DF( 1, 1 ) = ( ( tmp21 ) ? ( tmp70 ) : ( -tmp35 * tmp67 + tmp36 + tmp54 * tmp73 ) );
      DF( 1, 2 ) = ( ( tmp21 ) ? ( 0 ) : ( tmp54 * tmp77 ) );
      DF( 2, 0 ) = ( ( tmp21 ) ? ( 0 ) : ( -tmp48 * tmp78 + tmp57 * tmp79 - tmp61 * tmp80 ) );
      DF( 2, 1 ) = ( ( tmp21 ) ? ( 0 ) : ( -tmp69 * tmp78 + tmp71 * tmp79 - tmp72 * tmp80 ) );
      DF( 2, 2 ) = ( ( tmp21 ) ? ( 1 ) : ( -tmp34 * tmp74 + tmp51 * tmp76 - tmp78 * xold_2 ) );

      WALBERLA_ASSERT( !std::isnan( DF( 0, 0 ) ), "The torus map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 0, 1 ) ), "The torus map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 0, 2 ) ), "The torus map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 1, 0 ) ), "The torus map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 1, 1 ) ), "The torus map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 1, 2 ) ), "The torus map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 2, 0 ) ), "The torus map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 2, 1 ) ), "The torus map produces NaNs :( again." )
      WALBERLA_ASSERT( !std::isnan( DF( 2, 2 ) ), "The torus map produces NaNs :( again." )

      return DF.determinant();
   }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override final
   {
      WALBERLA_ABORT( "Serialization not implemented for TorusMap" );
   }

   /// \brief Applies the torus map to the SetupPrimitiveStorage.
   ///
   /// As unstructured base mesh, please use the mesh generator MeshInfo::meshTorus() with the same parameters for the torus.
   ///
   /// \param setupStorage the SetupPrimitiveStorage instance
   /// \param toroidalResolution number of prisms in toroidal direction (along the ring) in a complete (360 degree) ring
   /// \param poloidalResolution number of vertices on the boundary of a slice through the tube
   /// \param radiusOriginToCenterOfTube distance from origin to the center of the tube
   /// \param tubeLayerRadii list of radii of layers of the sliced tube - the last element defines the actual radius of the tube
   /// \param toroidalStartAngle angle (in radians) by which the domain shall be rotated about the z-axis
   /// \param poloidalStartAngle angle (in radians) by which the domain shall be rotated about the ring through the center of the tube
   ///
   static void setMap( SetupPrimitiveStorage& setupStorage,
                       uint_t                 toroidalResolution,
                       uint_t                 poloidalResolution,
                       real_t                 radiusOriginToCenterOfTube,
                       std::vector< real_t >  tubeLayerRadii,
                       real_t                 toroidalStartAngle,
                       real_t                 poloidalStartAngle )
   {
      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(),
                                      std::make_shared< TorusMap >( cell,
                                                                    setupStorage,
                                                                    toroidalResolution,
                                                                    poloidalResolution,
                                                                    radiusOriginToCenterOfTube,
                                                                    tubeLayerRadii,
                                                                    toroidalStartAngle,
                                                                    poloidalStartAngle ) );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(),
                                      std::make_shared< TorusMap >( face,
                                                                    setupStorage,
                                                                    toroidalResolution,
                                                                    poloidalResolution,
                                                                    radiusOriginToCenterOfTube,
                                                                    tubeLayerRadii,
                                                                    toroidalStartAngle,
                                                                    poloidalStartAngle ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(),
                                      std::make_shared< TorusMap >( edge,
                                                                    setupStorage,
                                                                    toroidalResolution,
                                                                    poloidalResolution,
                                                                    radiusOriginToCenterOfTube,
                                                                    tubeLayerRadii,
                                                                    toroidalStartAngle,
                                                                    poloidalStartAngle ) );
      }

      for ( auto it : setupStorage.getVertices() )
      {
         Vertex& vertex = *it.second;
         setupStorage.setGeometryMap( vertex.getID(),
                                      std::make_shared< TorusMap >( vertex,
                                                                    setupStorage,
                                                                    toroidalResolution,
                                                                    poloidalResolution,
                                                                    radiusOriginToCenterOfTube,
                                                                    tubeLayerRadii,
                                                                    toroidalStartAngle,
                                                                    poloidalStartAngle ) );
      }
   }

   /** @name 2D methods
  *    methods for 2D (class only provides a pseudo-implementation to satisfy requirements of base class)
  */
   ///@{
   void evalDF( const Point3D& x, Matrix2r& DFx ) const override final
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFx );
      WALBERLA_ABORT( "TorusMap::evalDF unimplemented for 2D!" );
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const override final
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFinvx );
      WALBERLA_ABORT( "TorusMap::evalDFinv unimplemented for 2D!" );
   }
   ///@}

 private:
   real_t angle( Point3D a, Point3D b ) const { return std::acos( a.dot( b ) / ( a.norm() * b.norm() ) ); }

   void identifyPrism( const Cell& cell )
   {
      tubeLayerRadiiBack_ = tubeLayerRadii_.back();

      toroidalAngleIncrement_ = 2 * pi / real_c( toroidalResolution );
      poloidalAngleIncrement_ = 2 * pi / real_c( poloidalResolution );

      auto    coords = cell.getCoordinates();
      Point3D centroid( 0, 0, 0 );
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

      auto    C                     = torusCoordinates( radiusOriginToCenterOfTube_, 0, toroidalAngle, 0 );
      Point3D centroidTrafoToOrigin = centroid - C;
      centroidTrafoToOrigin         = Point3D(
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

   uint_t                toroidalResolution;
   uint_t                poloidalResolution;
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

   real_t tubeLayerRadiiBack_;
};

} // end of namespace hyteg
