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

/// The class implements a mapping of the torus mesh to a tokamak domain.
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
               real_t                       r1,
               real_t                       r2 )
   : numToroidalSlices_( numToroidalSlices )
   , numPoloidalSlices_( numPoloidalSlices )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
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
               uint_t                       numToroidalSlices,
               uint_t                       numPoloidalSlices,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r1,
               real_t                       r2 )
   : numToroidalSlices_( numToroidalSlices )
   , numPoloidalSlices_( numPoloidalSlices )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
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
               uint_t                       numToroidalSlices,
               uint_t                       numPoloidalSlices,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r1,
               real_t                       r2 )
   : numToroidalSlices_( numToroidalSlices )
   , numPoloidalSlices_( numPoloidalSlices )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
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
               uint_t                       numToroidalSlices,
               uint_t                       numPoloidalSlices,
               real_t                       radiusOriginToCenterOfTube,
               std::vector< real_t >        tubeLayerRadii,
               real_t                       toroidalStartAngle,
               real_t                       poloidalStartAngle,
               real_t                       delta,
               real_t                       r1,
               real_t                       r2 )
   : numToroidalSlices_( numToroidalSlices )
   , numPoloidalSlices_( numPoloidalSlices )
   , radiusOriginToCenterOfTube_( radiusOriginToCenterOfTube )
   , tubeLayerRadii_( tubeLayerRadii )
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

   TokamakMap( walberla::mpi::RecvBuffer& recvBuffer ) { WALBERLA_ABORT( "Deserialization not implemented for TokamakMap" ); }

   void evalF( const Point3D& xold, Point3D& xnew ) const override final
   {
      // generated with data/codegen/geometry/TokamakMap.py

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
      auto tmp9  = -poloidalStartAngle_ + tmp8;
      auto tmp10 = ( ( poloidalStartAngle_ - tmp8 > 0 ) ? ( tmp9 + 6.2831853071795862 ) : ( tmp9 ) );
      auto tmp11 = poloidalStartAngle_ + tmp10;
      auto tmp12 = sin( tmp11 );
      auto tmp13 = 0.5 * poloidalAngleIncrement_;
      auto tmp14 = poloidalAngleIncrement_ * poloidalPrism_;
      auto tmp15 = sqrt( pow( tmp4, 2 ) + pow( tmp6, 2 ) + pow( xold_2, 2 ) ) *
                   sin( -tmp10 + tmp13 + ( ( tmp14 < 0 ) ? ( tmp14 + 6.2831853071795862 ) : ( tmp14 ) ) + 1.5707963267948966 ) /
                   ( tubeLayerRadiiBack_ * sin( tmp13 - 1.5707963267948966 ) );
      auto tmp16 = -r1_ * tmp15 * cos( tmp11 + tmp12 * asin( delta_ ) ) + radiusOriginToCenterOfTube_;
      xnew[0]    = tmp16 * tmp5;
      xnew[1]    = tmp16 * tmp7;
      xnew[2]    = -r2_ * tmp12 * tmp15;
   }

   real_t evalDF( const Point3D& xold, Matrix3r& DF ) const override final
   {
      // generated with data/codegen/geometry/TokamakMap.py

      auto xold_0 = xold[0];
      auto xold_1 = xold[1];
      auto xold_2 = xold[2];
      auto tmp0   = pow( xold_0, 2 );
      auto tmp1   = pow( xold_1, 2 );
      auto tmp2   = tmp0 + tmp1;
      auto tmp3   = pow( tmp2, -1.0 / 2.0 );
      auto tmp4   = asin( delta_ );
      auto tmp5   = radiusOriginToCenterOfTube_ * tmp3;
      auto tmp6   = 0.5 * toroidalAngleIncrement_;
      auto tmp7   = 1.0 / sin( tmp6 - 1.5707963267948966 );
      auto tmp8 =
          tmp6 + toroidalAngleIncrement_ * toroidalPrism_ + toroidalStartAngle_ - atan2( xold_1, xold_0 ) + 1.5707963267948966;
      auto tmp9  = tmp7 * sin( tmp8 );
      auto tmp10 = -tmp5 * xold_0 - tmp9 * xold_0;
      auto tmp11 = tmp10 * tmp3;
      auto tmp12 = -tmp5 * xold_1 - tmp9 * xold_1;
      auto tmp13 = tmp12 * tmp3;
      auto tmp14 = tmp11 * xold_0 + tmp13 * xold_1;
      auto tmp15 = atan2( xold_2, tmp14 );
      auto tmp16 = -poloidalStartAngle_ + tmp15;
      auto tmp17 = ( ( poloidalStartAngle_ - tmp15 > 0 ) ? ( tmp16 + 6.2831853071795862 ) : ( tmp16 ) );
      auto tmp18 = poloidalStartAngle_ + tmp17;
      auto tmp19 = sin( tmp18 );
      auto tmp20 = tmp18 + tmp19 * tmp4;
      auto tmp21 = cos( tmp20 );
      auto tmp22 = pow( xold_2, 2 );
      auto tmp23 = sqrt( pow( tmp10, 2 ) + pow( tmp12, 2 ) + tmp22 );
      auto tmp24 = 1.0 / tubeLayerRadiiBack_;
      auto tmp25 = 0.5 * poloidalAngleIncrement_;
      auto tmp26 = 1.0 / sin( tmp25 - 1.5707963267948966 );
      auto tmp27 = poloidalAngleIncrement_ * poloidalPrism_;
      auto tmp28 = -tmp17 + tmp25 + ( ( tmp27 < 0 ) ? ( tmp27 + 6.2831853071795862 ) : ( tmp27 ) ) + 1.5707963267948966;
      auto tmp29 = tmp24 * tmp26 * sin( tmp28 );
      auto tmp30 = tmp23 * tmp29;
      auto tmp31 = r1_ * tmp30;
      auto tmp32 = radiusOriginToCenterOfTube_ - tmp21 * tmp31;
      auto tmp33 = tmp3 * tmp32;
      auto tmp34 = pow( tmp2, -3.0 / 2.0 );
      auto tmp35 = tmp0 * tmp34;
      auto tmp36 = r1_ * tmp21;
      auto tmp37 = xold_0 * xold_1;
      auto tmp38 = tmp34 * tmp37;
      auto tmp39 = radiusOriginToCenterOfTube_ * tmp38;
      auto tmp40 = 2 * tmp39;
      auto tmp41 = tmp7 * cos( tmp8 ) / tmp2;
      auto tmp42 = tmp1 * tmp41;
      auto tmp43 = ( 1.0 / 2.0 ) * tmp12;
      auto tmp44 = radiusOriginToCenterOfTube_ * tmp35;
      auto tmp45 = tmp37 * tmp41;
      auto tmp46 = 2 * tmp45;
      auto tmp47 = -2 * tmp5 - 2 * tmp9;
      auto tmp48 = ( 1.0 / 2.0 ) * tmp10;
      auto tmp49 = tmp29 / tmp23;
      auto tmp50 = tmp49 * ( tmp43 * ( tmp40 - 2 * tmp42 ) + tmp48 * ( 2 * tmp44 - tmp46 + tmp47 ) );
      auto tmp51 = tmp3 * xold_1;
      auto tmp52 = -tmp5 - tmp9;
      auto tmp53 = tmp3 * xold_0;
      auto tmp54 = 1.0 / ( pow( tmp14, 2 ) + tmp22 );
      auto tmp55 = tmp54 * xold_2;
      auto tmp56 =
          tmp55 * ( -tmp10 * tmp35 + tmp11 - tmp12 * tmp38 + tmp51 * ( tmp39 - tmp42 ) + tmp53 * ( tmp44 - tmp45 + tmp52 ) );
      auto tmp57 = tmp23 * tmp24 * tmp26 * cos( tmp28 );
      auto tmp58 = tmp36 * tmp57;
      auto tmp59 = cos( tmp18 );
      auto tmp60 = tmp4 * tmp59;
      auto tmp61 = tmp31 * sin( tmp20 );
      auto tmp62 = -tmp36 * tmp50 - tmp56 * tmp58 + tmp61 * ( -tmp56 * tmp60 - tmp56 );
      auto tmp63 = -tmp32 * tmp38;
      auto tmp64 = tmp0 * tmp41;
      auto tmp65 = tmp1 * tmp34;
      auto tmp66 = radiusOriginToCenterOfTube_ * tmp65;
      auto tmp67 = tmp43 * ( tmp46 + tmp47 + 2 * tmp66 ) + tmp48 * ( tmp40 + 2 * tmp64 );
      auto tmp68 = tmp36 * tmp49;
      auto tmp69 =
          tmp55 * ( -tmp10 * tmp38 - tmp12 * tmp65 + tmp13 + tmp51 * ( tmp45 + tmp52 + tmp66 ) + tmp53 * ( tmp39 + tmp64 ) );
      auto tmp70 = -tmp58 * tmp69 + tmp61 * ( -tmp60 * tmp69 - tmp69 ) - tmp67 * tmp68;
      auto tmp71 = tmp14 * tmp54;
      auto tmp72 = tmp57 * tmp71;
      auto tmp73 = tmp36 * tmp72 + tmp61 * ( tmp60 * tmp71 + tmp71 ) - tmp68 * xold_2;
      auto tmp74 = r2_ * tmp19;
      auto tmp75 = r2_ * tmp30 * tmp59;
      auto tmp76 = tmp57 * tmp74;
      auto tmp77 = tmp49 * tmp74;
      DF( 0, 0 ) = -tmp32 * tmp35 + tmp33 + tmp53 * tmp62;
      DF( 0, 1 ) = tmp53 * tmp70 + tmp63;
      DF( 0, 2 ) = tmp53 * tmp73;
      DF( 1, 0 ) = tmp51 * tmp62 + tmp63;
      DF( 1, 1 ) = -tmp32 * tmp65 + tmp33 + tmp51 * tmp70;
      DF( 1, 2 ) = tmp51 * tmp73;
      DF( 2, 0 ) = -tmp50 * tmp74 + tmp56 * tmp75 - tmp56 * tmp76;
      DF( 2, 1 ) = -tmp67 * tmp77 + tmp69 * tmp75 - tmp69 * tmp76;
      DF( 2, 2 ) = -tmp71 * tmp75 + tmp72 * tmp74 - tmp77 * xold_2;
      return DF.determinant();
   }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override final
   {
      WALBERLA_ABORT( "Serialization not implemented for TokamakMap" );
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
   /// \param setupStorage the SetupPrimitiveStorage instance
   /// \param numToroidalSlices number of prisms in toroidal direction (along the ring)
   /// \param numPoloidalSlices number of vertices on the boundary of a slice through the tube
   /// \param radiusOriginToCenterOfTube distance from origin to the center of the tube
   /// \param tubeLayerRadii list of radii of layers of the sliced tube - the last element defines the actual radius of the tube
   /// \param toroidalStartAngle angle (in radians) by which the domain shall be rotated about the z-axis
   /// \param poloidalStartAngle angle (in radians) by which the domain shall be rotated about the ring through the center of the tube
   /// \param delta triangularity parameter of the tokamak
   /// \param r1 semi-minor axis radius
   /// \param r2 semi-major axis radius
   ///
   static void setMap( SetupPrimitiveStorage& setupStorage,
                       uint_t                 numToroidalSlices,
                       uint_t                 numPoloidalSlices,
                       real_t                 radiusOriginToCenterOfTube,
                       std::vector< real_t >  tubeLayerRadii,
                       real_t                 toroidalStartAngle,
                       real_t                 poloidalStartAngle,
                       real_t                 delta,
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
                                                                      r1,
                                                                      r2 ) );
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
      WALBERLA_ABORT( "TokamakMap::evalDF unimplemented for 2D!" );
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const override final
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
      Point3D centroidTrafoToOrigin = centroid - C;
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
   real_t r1_;
   real_t r2_;

   real_t tubeLayerRadiiBack_;
};

} // end of namespace hyteg
