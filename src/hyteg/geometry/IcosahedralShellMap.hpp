/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr, Andreas Burkhart.
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

#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/types/Matrix.hpp"

#include "GeometryMap.hpp"

// #define SHELL_MAP_LOG( STR ) WALBERLA_LOG_INFO_ON_ROOT( STR );
#define SHELL_MAP_LOG( STR )

namespace hyteg {

using walberla::real_c;

/// Class providing geometry mapping for a facetted isosahedral shell
///
/// This geometry map provides a blending operation for a base mesh generated
/// with the inline meshSphericalShell generator. Geometric nodes on refined
/// hierarchy levels are projected onto spherical layers.
///
/// Note that although the mesh nodes are arranged on radial shells, the
/// nodes are not aligned on radial beams through the origin.
/// The more expensive IcosahedralShellAlignedMap does have that feature.
/// (For a visual comparison in 2D have a look at the documentation of the AnnulusMap
/// and the AnnulusAlignedMap.)
///
class IcosahedralShellMap : public GeometryMap
{
 public:
   IcosahedralShellMap( const Cell& cell, const SetupPrimitiveStorage& storage )
   {
      SHELL_MAP_LOG( "---------------------------------------------" );
      SHELL_MAP_LOG( "Initialising Shell map for cellID: " << cell.getID() );
      SHELL_MAP_LOG( "---------------------------------------------" );
      classifyVertices( cell, storage );
   }

   IcosahedralShellMap( const Face& face, const SetupPrimitiveStorage& storage )
   {
      SHELL_MAP_LOG( "---------------------------------------------" );
      SHELL_MAP_LOG( "Initialising Shell map for faceID: " << face.getID() );
      SHELL_MAP_LOG( "---------------------------------------------" );

      std::vector< PrimitiveID > neighborCells;
      face.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *storage.getCell( neighborCells[0] );
      classifyVertices( cell, storage );
   }

   IcosahedralShellMap( const Edge& edge, const SetupPrimitiveStorage& storage )
   {
      SHELL_MAP_LOG( "---------------------------------------------" );
      SHELL_MAP_LOG( "Initialising Shell map for edgeID: " << edge.getID() );
      SHELL_MAP_LOG( "---------------------------------------------" );

      std::vector< PrimitiveID > neighborCells;
      edge.getNeighborCells( neighborCells );
      WALBERLA_ASSERT_GREATER( neighborCells.size(), 0 );
      const Cell& cell = *storage.getCell( neighborCells[0] );
      classifyVertices( cell, storage );
   }

   IcosahedralShellMap( walberla::mpi::RecvBuffer& recvBuffer )
   {
      recvBuffer >> rayVertex_[0];
      recvBuffer >> rayVertex_[1];
      recvBuffer >> rayVertex_[2];

      recvBuffer >> refVertex_[0];
      recvBuffer >> refVertex_[1];
      recvBuffer >> refVertex_[2];

      recvBuffer >> thrVertex_[0];
      recvBuffer >> thrVertex_[1];
      recvBuffer >> thrVertex_[2];

      recvBuffer >> forVertex_[0];
      recvBuffer >> forVertex_[1];
      recvBuffer >> forVertex_[2];

      recvBuffer >> radRefVertex_;
      recvBuffer >> radRayVertex_;

      recvBuffer >> prismNormal_;
   }

   void evalF( const Point3D& xold, Point3D& xnew ) const override
   {
      // determine barycentric coordinate w.r.t. vertex refVertex_
      real_t tmp0  = -rayVertex_[2];
      real_t tmp1  = refVertex_[2] + tmp0;
      real_t tmp2  = -rayVertex_[0];
      real_t tmp3  = thrVertex_[0] + tmp2;
      real_t tmp4  = -rayVertex_[1];
      real_t tmp5  = forVertex_[1] + tmp4;
      real_t tmp6  = tmp3 * tmp5;
      real_t tmp7  = refVertex_[1] + tmp4;
      real_t tmp8  = thrVertex_[2] + tmp0;
      real_t tmp9  = forVertex_[0] + tmp2;
      real_t tmp10 = tmp8 * tmp9;
      real_t tmp11 = refVertex_[0] + tmp2;
      real_t tmp12 = thrVertex_[1] + tmp4;
      real_t tmp13 = forVertex_[2] + tmp0;
      real_t tmp14 = tmp12 * tmp13;
      real_t tmp15 = tmp13 * tmp3;
      real_t tmp16 = tmp12 * tmp9;
      real_t tmp17 = tmp5 * tmp8;
      real_t tmp18 = tmp0 + xold[2];
      real_t tmp19 = tmp4 + xold[1];
      real_t tmp20 = tmp2 + xold[0];

      real_t volT = -tmp1 * tmp16 + tmp1 * tmp6 + tmp10 * tmp7 + tmp11 * tmp14 - tmp11 * tmp17 - tmp15 * tmp7;
      real_t volX = tmp10 * tmp19 + tmp14 * tmp20 - tmp15 * tmp19 - tmp16 * tmp18 - tmp17 * tmp20 + tmp18 * tmp6;
      real_t bary = std::abs( volX / volT );

      // compute new coordinates
      real_t oldRad = std::sqrt( xold[0] * xold[0] + xold[1] * xold[1] + xold[2] * xold[2] );
      real_t newRad = radRayVertex_ + bary * ( radRefVertex_ - radRayVertex_ );
      xnew[0]       = xold[0] * newRad / oldRad;
      xnew[1]       = xold[1] * newRad / oldRad;
      xnew[2]       = xold[2] * newRad / oldRad;

      // SHELL_MAP_LOG( "Mapped: " << xold << " --> " << xnew );
   }

   void evalFinv( const Point3D& xPhys, Point3D& xComp ) const override
   {
      // calculating the intersection point of the prism-parallel plane that contains xComp
      // with the line spanned by mu * xPhys

      // plane through xComp given by all points p with ( p - A' ) * n = 0
      // where A' is the point on a ray that has distance norm(xPhys) from origin
      // (lies on the sphere of xPhys)
      // therefore: A' = A * (norm(xPhys) / norm(A)) where A is any point on a ray (we use A = refVertex_)

      // since xComp = mu * xPhys, we insert in plane eq. and solve for mu to find xComp
      // solve for mu:
      //    ((mu * xPhys) - (A * (norm(xPhys) / norm(A)))) * n = 0
      // => mu = (norm(xPhys) / norm(A)) * ( (A * n) / (xPhys * n) )
      // let's do that

      const auto A  = refVertex_;
      const auto mu = ( xPhys.norm() / A.norm() ) * ( A.dot( prismNormal_ ) / xPhys.dot( prismNormal_ ) );
      xComp         = mu * xPhys;
   }

   void evalDF( const Point3D& x, Matrix2r& DFx ) const final override
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFx );
      WALBERLA_ABORT( "IcosahedralMap::evalDF unimplemented for 2D!" );
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const final override
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFinvx );
      WALBERLA_ABORT( "IcosahedralMap::evalDFinv unimplemented for 2D!" );
   }

   real_t evalDF( const Point3D& x, Matrix3r& DFx ) const final override
   {
      // real_t tmp0 = pow(x[0], 2);
      real_t tmp0  = x[0] * x[0];
      real_t tmp1  = rayVertex_[2] - refVertex_[2];
      real_t tmp2  = rayVertex_[0] - thrVertex_[0];
      real_t tmp3  = rayVertex_[1] - forVertex_[1];
      real_t tmp4  = tmp2 * tmp3;
      real_t tmp5  = rayVertex_[1] - refVertex_[1];
      real_t tmp6  = rayVertex_[0] - forVertex_[0];
      real_t tmp7  = rayVertex_[2] - thrVertex_[2];
      real_t tmp8  = tmp6 * tmp7;
      real_t tmp9  = rayVertex_[0] - refVertex_[0];
      real_t tmp10 = rayVertex_[1] - thrVertex_[1];
      real_t tmp11 = rayVertex_[2] - forVertex_[2];
      real_t tmp12 = tmp10 * tmp11;
      real_t tmp13 = tmp11 * tmp2;
      real_t tmp14 = tmp10 * tmp6;
      real_t tmp15 = tmp3 * tmp7;
      real_t tmp16 = -tmp1 * tmp14 + tmp1 * tmp4 + tmp12 * tmp9 - tmp13 * tmp5 - tmp15 * tmp9 + tmp5 * tmp8;
      real_t tmp17 = radRayVertex_ - radRefVertex_;
      real_t tmp18 = rayVertex_[2] - x[2];
      real_t tmp19 = rayVertex_[1] - x[1];
      real_t tmp20 = rayVertex_[0] - x[0];
      real_t tmp21 = radRayVertex_ * tmp16 -
                     tmp17 * ( tmp12 * tmp20 - tmp13 * tmp19 - tmp14 * tmp18 - tmp15 * tmp20 + tmp18 * tmp4 + tmp19 * tmp8 );
      // real_t tmp22 = pow(x[1], 2);
      // real_t tmp23 = pow(x[2], 2);
      real_t tmp22 = x[1] * x[1];
      real_t tmp23 = x[2] * x[2];
      real_t tmp24 = tmp0 + tmp22 + tmp23;
      real_t tmp25 = tmp17 * ( tmp12 - tmp15 );
      real_t tmp26 = real_c( 1.0 / ( tmp16 * pow( tmp24, 3.0 / 2.0 ) ) );
      real_t tmp27 = tmp13 - tmp8;
      real_t tmp28 = tmp17 * tmp24;
      real_t tmp29 = tmp21 * x[1] + tmp27 * tmp28;
      real_t tmp30 = tmp26 * x[0];
      real_t tmp31 = -tmp14 + tmp4;
      real_t tmp32 = -tmp21 * x[2] + tmp28 * tmp31;
      real_t tmp33 = -tmp21 * x[0] + tmp24 * tmp25;
      real_t tmp34 = tmp26 * x[1];
      real_t tmp35 = tmp26 * x[2];

      DFx( 0, 0 ) = tmp26 * ( -tmp0 * tmp21 + tmp24 * ( tmp21 + tmp25 * x[0] ) );
      DFx( 0, 1 ) = -tmp29 * tmp30;
      DFx( 0, 2 ) = tmp30 * tmp32;
      DFx( 1, 0 ) = tmp33 * tmp34;
      DFx( 1, 1 ) = tmp26 * ( -tmp21 * tmp22 + tmp24 * ( -tmp17 * tmp27 * x[1] + tmp21 ) );
      DFx( 1, 2 ) = tmp32 * tmp34;
      DFx( 2, 0 ) = tmp33 * tmp35;
      DFx( 2, 1 ) = -tmp29 * tmp35;
      DFx( 2, 2 ) = tmp26 * ( -tmp21 * tmp23 + tmp24 * ( tmp17 * tmp31 * x[2] + tmp21 ) );

      return DFx( 0, 0 ) * DFx( 1, 1 ) * DFx( 2, 2 ) - DFx( 0, 0 ) * DFx( 2, 1 ) * DFx( 1, 2 ) -
             DFx( 1, 0 ) * DFx( 0, 1 ) * DFx( 2, 2 ) + DFx( 1, 0 ) * DFx( 2, 1 ) * DFx( 0, 2 ) +
             DFx( 2, 0 ) * DFx( 0, 1 ) * DFx( 1, 2 ) - DFx( 2, 0 ) * DFx( 1, 1 ) * DFx( 0, 2 );
   };

   void evalDFinv( const Point3D& x, Matrix3r& DFx ) const final override
   {
      Matrix3r tmp;
      evalDF( x, tmp );

      DFx = tmp.inverse();
   }

   void evalDFinvDF( const Point3D& x, Matrixr< 3, 9 >& DFinvDFx ) const override final
   {
      const real_t tmp0 = ( x[0] * x[0] );
      const real_t tmp1 = ( x[1] * x[1] );
      const real_t tmp2 = ( x[2] * x[2] );
      const real_t tmp3 = tmp0 + tmp1 + tmp2;
      const real_t tmp4 = std::pow( tmp3, -1.0 / 2.0 );
      const real_t tmp5 = std::pow(
          ( rayVertex_[0] * rayVertex_[0] ) + ( rayVertex_[1] * rayVertex_[1] ) + ( rayVertex_[2] * rayVertex_[2] ), 1.0 / 2.0 );
      const real_t tmp6 = -tmp5 + std::pow( ( refVertex_[0] * refVertex_[0] ) + ( refVertex_[1] * refVertex_[1] ) +
                                                ( refVertex_[2] * refVertex_[2] ),
                                            1.0 / 2.0 );
      const real_t tmp7 = forVertex_[0] * rayVertex_[1];
      const real_t tmp8 = forVertex_[1] * thrVertex_[0];
      const real_t tmp9 = rayVertex_[0] * thrVertex_[1];
      const real_t tmp10 =
          -tmp7 - tmp8 - tmp9 + forVertex_[0] * thrVertex_[1] + forVertex_[1] * rayVertex_[0] + rayVertex_[1] * thrVertex_[0];
      const real_t tmp11 = forVertex_[0] * rayVertex_[2];
      const real_t tmp12 = tmp11 * thrVertex_[1];
      const real_t tmp13 = forVertex_[0] * thrVertex_[2];
      const real_t tmp14 = rayVertex_[0] * thrVertex_[2];
      const real_t tmp15 = tmp14 * forVertex_[1];
      const real_t tmp16 = forVertex_[1] * rayVertex_[2];
      const real_t tmp17 = forVertex_[2] * rayVertex_[0];
      const real_t tmp18 = forVertex_[2] * thrVertex_[0];
      const real_t tmp19 = tmp18 * rayVertex_[1];
      const real_t tmp20 = forVertex_[2] * thrVertex_[1];
      const real_t tmp21 = rayVertex_[1] * thrVertex_[2];
      const real_t tmp22 = rayVertex_[2] * thrVertex_[0];
      const real_t tmp23 =
          1.0 /
          ( -tmp12 - tmp13 * refVertex_[1] - tmp15 - tmp16 * refVertex_[0] - tmp17 * refVertex_[1] - tmp19 -
            tmp20 * refVertex_[0] - tmp21 * refVertex_[0] - tmp22 * refVertex_[1] - tmp7 * refVertex_[2] - tmp8 * refVertex_[2] -
            tmp9 * refVertex_[2] + forVertex_[0] * rayVertex_[1] * thrVertex_[2] + forVertex_[0] * rayVertex_[2] * refVertex_[1] +
            forVertex_[0] * refVertex_[2] * thrVertex_[1] + forVertex_[1] * rayVertex_[0] * refVertex_[2] +
            forVertex_[1] * rayVertex_[2] * thrVertex_[0] + forVertex_[1] * refVertex_[0] * thrVertex_[2] +
            forVertex_[2] * rayVertex_[0] * thrVertex_[1] + forVertex_[2] * rayVertex_[1] * refVertex_[0] +
            forVertex_[2] * refVertex_[1] * thrVertex_[0] + rayVertex_[0] * refVertex_[1] * thrVertex_[2] +
            rayVertex_[1] * refVertex_[2] * thrVertex_[0] + rayVertex_[2] * refVertex_[0] * thrVertex_[1] );
      const real_t tmp24 = tmp23 * tmp6;
      const real_t tmp25 = tmp24 * ( tmp11 * x[1] - tmp12 + tmp13 * rayVertex_[1] - tmp13 * x[1] + tmp14 * x[1] - tmp15 -
                                     tmp16 * x[0] + tmp17 * thrVertex_[1] - tmp17 * x[1] + tmp18 * x[1] - tmp19 - tmp20 * x[0] -
                                     tmp21 * x[0] + tmp22 * forVertex_[1] - tmp22 * x[1] - tmp7 * x[2] - tmp8 * x[2] -
                                     tmp9 * x[2] + forVertex_[0] * thrVertex_[1] * x[2] + forVertex_[1] * rayVertex_[0] * x[2] +
                                     forVertex_[1] * thrVertex_[2] * x[0] + forVertex_[2] * rayVertex_[1] * x[0] +
                                     rayVertex_[1] * thrVertex_[0] * x[2] + rayVertex_[2] * thrVertex_[1] * x[0] ) +
                           tmp5;
      const real_t tmp26 = std::pow( tmp3, -3.0 / 2.0 );
      const real_t tmp27 = tmp25 * tmp26;
      const real_t tmp28 = tmp27 * x[1];
      const real_t tmp29 = tmp28 * x[2];
      const real_t tmp30 = -tmp10 * tmp23 * tmp4 * tmp6 * x[1] + tmp29;
      const real_t tmp31 = -tmp30;
      const real_t tmp32 = tmp11 - tmp13 + tmp14 - tmp17 + tmp18 - tmp22;
      const real_t tmp33 = tmp24 * tmp26;
      const real_t tmp34 = tmp32 * tmp33;
      const real_t tmp35 = tmp34 * x[0];
      const real_t tmp36 = tmp35 * x[2];
      const real_t tmp37 = std::pow( tmp3, -5.0 / 2.0 );
      const real_t tmp38 = 3 * tmp25 * tmp37;
      const real_t tmp39 = tmp38 * x[0];
      const real_t tmp40 = x[1] * x[2];
      const real_t tmp41 = -tmp39 * tmp40;
      const real_t tmp42 =
          -tmp16 - tmp20 - tmp21 + forVertex_[1] * thrVertex_[2] + forVertex_[2] * rayVertex_[1] + rayVertex_[2] * thrVertex_[1];
      const real_t tmp43 = tmp33 * tmp42;
      const real_t tmp44 = tmp40 * tmp43 + tmp41;
      const real_t tmp45 = tmp36 + tmp44;
      const real_t tmp46 = tmp31 * tmp45;
      const real_t tmp47 = -tmp23 * tmp32 * tmp4 * tmp6 * x[2] + tmp29;
      const real_t tmp48 = tmp10 * tmp33;
      const real_t tmp49 = tmp48 * x[0];
      const real_t tmp50 = tmp49 * x[1];
      const real_t tmp51 = tmp44 + tmp50;
      const real_t tmp52 = -tmp51;
      const real_t tmp53 = tmp47 * tmp52;
      const real_t tmp54 = tmp24 * tmp4;
      const real_t tmp55 = tmp32 * tmp54;
      const real_t tmp56 = tmp25 * tmp4;
      const real_t tmp57 = -tmp1 * tmp27 + tmp55 * x[1] + tmp56;
      const real_t tmp58 = -tmp23 * tmp4 * tmp42 * tmp6;
      const real_t tmp59 = tmp49 * x[2];
      const real_t tmp60 = tmp27 * x[0];
      const real_t tmp61 = -3 * tmp2 * tmp25 * tmp37 * x[0] + tmp60;
      const real_t tmp62 = tmp2 * tmp43 + tmp58 + tmp59 + tmp61;
      const real_t tmp63 = -tmp62;
      const real_t tmp64 = tmp57 * tmp63;
      const real_t tmp65 = tmp10 * tmp54 * x[2] - tmp2 * tmp27 + tmp56;
      const real_t tmp66 = tmp35 * x[1];
      const real_t tmp67 = -tmp1 * tmp39 + tmp60;
      const real_t tmp68 = tmp1 * tmp43 + tmp58 + tmp66 + tmp67;
      const real_t tmp69 = -tmp68;
      const real_t tmp70 = tmp65 * tmp69;
      const real_t tmp71 = tmp60 * x[1];
      const real_t tmp72 = -tmp23 * tmp32 * tmp4 * tmp6 * x[0] + tmp71;
      const real_t tmp73 = -tmp72;
      const real_t tmp74 = tmp60 * x[2];
      const real_t tmp75 = -tmp23 * tmp4 * tmp42 * tmp6 * x[2] + tmp74;
      const real_t tmp76 = -tmp75;
      const real_t tmp77 = -tmp23 * tmp4 * tmp42 * tmp6 * x[1] + tmp71;
      const real_t tmp78 = -tmp77;
      const real_t tmp79 = -tmp10 * tmp23 * tmp4 * tmp6 * x[0] + tmp74;
      const real_t tmp80 = -tmp79;
      const real_t tmp81 = -tmp47;
      const real_t tmp82 = tmp65 * tmp73;
      const real_t tmp83 = tmp57 * tmp80;
      const real_t tmp84 = -tmp0 * tmp27 + tmp42 * tmp54 * x[0] + tmp56;
      const real_t tmp85 = tmp31 * tmp81;
      const real_t tmp86 =
          tmp31 * tmp73 * tmp76 + tmp57 * tmp65 * tmp84 - tmp76 * tmp83 + tmp78 * tmp80 * tmp81 - tmp78 * tmp82 - tmp84 * tmp85;
      const real_t tmp87  = 1.0 / ( tmp86 );
      const real_t tmp88  = tmp57 * tmp65 - tmp85;
      const real_t tmp89  = 1.0 / ( tmp86 * tmp86 );
      const real_t tmp90  = tmp52 * tmp73;
      const real_t tmp91  = -tmp45;
      const real_t tmp92  = tmp80 * tmp91;
      const real_t tmp93  = 2 * tmp43;
      const real_t tmp94  = tmp93 * x[0];
      const real_t tmp95  = tmp27 * x[2];
      const real_t tmp96  = -3 * tmp0 * tmp25 * tmp37 * x[2] + tmp95;
      const real_t tmp97  = tmp94 * x[2] + tmp96;
      const real_t tmp98  = -tmp97;
      const real_t tmp99  = tmp31 * tmp73;
      const real_t tmp100 = -tmp0 * tmp38 * x[1] + tmp28;
      const real_t tmp101 = tmp100 + tmp94 * x[1];
      const real_t tmp102 = -tmp101;
      const real_t tmp103 = tmp80 * tmp81;
      const real_t tmp104 =
          -tmp0 * tmp93 + 2 * tmp23 * tmp4 * tmp42 * tmp6 + 3 * tmp25 * tmp37 * ( x[0] * x[0] * x[0] ) - 3 * tmp60;
      const real_t tmp105 = tmp104 * tmp47;
      const real_t tmp106 = tmp102 * tmp72;
      const real_t tmp107 = tmp79 * tmp98;
      const real_t tmp108 = tmp63 * tmp72;
      const real_t tmp109 = -tmp10 * tmp23 * tmp4 * tmp6;
      const real_t tmp110 = tmp0 * tmp33;
      const real_t tmp111 = tmp43 * x[0];
      const real_t tmp112 = tmp10 * tmp110 + tmp109 + tmp111 * x[2] + tmp96;
      const real_t tmp113 = -tmp112;
      const real_t tmp114 = tmp113 * tmp81;
      const real_t tmp115 = tmp69 * tmp79;
      const real_t tmp116 = -tmp55;
      const real_t tmp117 = tmp100 + tmp110 * tmp32 + tmp111 * x[1] + tmp116;
      const real_t tmp118 = -tmp117;
      const real_t tmp119 = tmp118 * tmp31;
      const real_t tmp120 = tmp117 * tmp65;
      const real_t tmp121 = tmp112 * tmp57;
      const real_t tmp122 = tmp57 * tmp65;
      const real_t tmp123 =
          tmp89 * ( -tmp102 * tmp103 - tmp104 * tmp122 - tmp105 * tmp31 - tmp106 * tmp65 - tmp107 * tmp57 - tmp108 * tmp78 -
                    tmp114 * tmp78 - tmp115 * tmp76 - tmp119 * tmp76 - tmp120 * tmp78 - tmp121 * tmp76 - tmp46 * tmp84 -
                    tmp53 * tmp84 - tmp64 * tmp84 - tmp70 * tmp84 - tmp76 * tmp90 - tmp78 * tmp92 - tmp98 * tmp99 );
      const real_t tmp124 = tmp80 * tmp81 - tmp82;
      const real_t tmp125 = -tmp83 + tmp99;
      const real_t tmp126 = tmp65 * tmp78;
      const real_t tmp127 = -tmp126 + tmp31 * tmp76;
      const real_t tmp128 = tmp65 * tmp84 - tmp76 * tmp80;
      const real_t tmp129 = tmp78 * tmp80;
      const real_t tmp130 = tmp31 * tmp84;
      const real_t tmp131 = tmp129 - tmp130;
      const real_t tmp132 = tmp78 * tmp81;
      const real_t tmp133 = tmp57 * tmp76;
      const real_t tmp134 = tmp132 - tmp133;
      const real_t tmp135 = tmp73 * tmp76;
      const real_t tmp136 = tmp135 - tmp81 * tmp84;
      const real_t tmp137 = tmp57 * tmp84 - tmp73 * tmp78;
      const real_t tmp138 = tmp34 * tmp40;
      const real_t tmp139 = -tmp1 * tmp38 * x[2] + tmp95;
      const real_t tmp140 = 2 * tmp138 + tmp139;
      const real_t tmp141 = tmp1 * tmp48 + tmp109 + tmp138 + tmp139;
      const real_t tmp142 = -tmp141;
      const real_t tmp143 = tmp142 * tmp47;
      const real_t tmp144 =
          -2 * tmp1 * tmp34 + 2 * tmp23 * tmp32 * tmp4 * tmp6 + 3 * tmp25 * tmp37 * ( x[1] * x[1] * x[1] ) - 3 * tmp28;
      const real_t tmp145 = tmp40 * tmp48;
      const real_t tmp146 = -3 * tmp2 * tmp25 * tmp37 * x[1] + tmp28;
      const real_t tmp147 = tmp116 + tmp145 + tmp146 + tmp2 * tmp34;
      const real_t tmp148 = -tmp147;
      const real_t tmp149 = tmp36 + tmp41 + tmp50;
      const real_t tmp150 = -tmp149;
      const real_t tmp151 = -tmp140;
      const real_t tmp152 = 2 * tmp66 + tmp67;
      const real_t tmp153 = -tmp152;
      const real_t tmp154 = tmp31 * tmp76;
      const real_t tmp155 = tmp144 * tmp79;
      const real_t tmp156 = tmp79 * tmp91;
      const real_t tmp157 = tmp148 * tmp72;
      const real_t tmp158 = tmp65 * tmp84;
      const real_t tmp159 = tmp57 * tmp84;
      const real_t tmp160 =
          tmp89 * ( -tmp103 * tmp69 - tmp118 * tmp122 - tmp119 * tmp47 - tmp126 * tmp152 - tmp129 * tmp151 - tmp130 * tmp140 -
                    tmp132 * tmp150 - tmp133 * tmp149 - tmp135 * tmp142 - tmp143 * tmp84 - tmp144 * tmp158 - tmp148 * tmp159 -
                    tmp153 * tmp154 - tmp155 * tmp76 - tmp156 * tmp57 - tmp157 * tmp78 - tmp70 * tmp72 - tmp91 * tmp99 );
      const real_t tmp161 = tmp142 * tmp76;
      const real_t tmp162 = 2 * tmp145 + tmp146;
      const real_t tmp163 = -tmp162;
      const real_t tmp164 = tmp163 * tmp47;
      const real_t tmp165 =
          2 * tmp10 * tmp23 * tmp4 * tmp6 - 2 * tmp2 * tmp48 + 3 * tmp25 * tmp37 * ( x[2] * x[2] * x[2] ) - 3 * tmp95;
      const real_t tmp166 = 2 * tmp59 + tmp61;
      const real_t tmp167 = -tmp166;
      const real_t tmp168 = tmp165 * tmp72;
      const real_t tmp169 = tmp52 * tmp72;
      const real_t tmp170 = tmp113 * tmp47;
      const real_t tmp171 =
          tmp89 * ( -tmp103 * tmp52 - tmp113 * tmp122 - tmp126 * tmp149 - tmp129 * tmp148 - tmp130 * tmp147 - tmp132 * tmp167 -
                    tmp133 * tmp166 - tmp135 * tmp163 - tmp142 * tmp158 - tmp150 * tmp154 - tmp159 * tmp165 - tmp161 * tmp79 -
                    tmp164 * tmp84 - tmp168 * tmp78 - tmp169 * tmp65 - tmp170 * tmp31 - tmp63 * tmp99 - tmp64 * tmp79 );
      DFinvDFx( 0, 0 ) = tmp123 * tmp88 + tmp87 * ( tmp46 + tmp53 + tmp64 + tmp70 );
      DFinvDFx( 0, 1 ) = tmp123 * tmp124 + tmp87 * ( tmp108 + tmp114 + tmp120 + tmp92 );
      DFinvDFx( 0, 2 ) = tmp123 * tmp125 + tmp87 * ( tmp115 + tmp119 + tmp121 + tmp90 );
      DFinvDFx( 1, 0 ) = tmp123 * tmp127 + tmp87 * ( tmp101 * tmp65 + tmp31 * tmp98 + tmp52 * tmp76 + tmp63 * tmp77 );
      DFinvDFx( 1, 1 ) = tmp123 * tmp128 + tmp87 * ( tmp104 * tmp65 + tmp107 + tmp112 * tmp76 + tmp63 * tmp84 );
      DFinvDFx( 1, 2 ) = tmp123 * tmp131 + tmp87 * ( tmp102 * tmp80 + tmp104 * tmp30 + tmp113 * tmp78 + tmp51 * tmp84 );
      DFinvDFx( 2, 0 ) = tmp123 * tmp134 + tmp87 * ( tmp102 * tmp81 + tmp57 * tmp97 + tmp69 * tmp75 + tmp78 * tmp91 );
      DFinvDFx( 2, 1 ) = tmp123 * tmp136 + tmp87 * ( tmp105 + tmp118 * tmp76 + tmp45 * tmp84 + tmp73 * tmp98 );
      DFinvDFx( 2, 2 ) = tmp123 * tmp137 + tmp87 * ( tmp104 * tmp57 + tmp106 + tmp117 * tmp78 + tmp69 * tmp84 );
      DFinvDFx( 0, 3 ) = tmp160 * tmp88 + tmp87 * ( tmp140 * tmp31 + tmp143 + tmp144 * tmp65 + tmp148 * tmp57 );
      DFinvDFx( 0, 4 ) = tmp124 * tmp160 + tmp87 * ( tmp150 * tmp81 + tmp151 * tmp80 + tmp152 * tmp65 + tmp157 );
      DFinvDFx( 0, 5 ) = tmp125 * tmp160 + tmp87 * ( tmp142 * tmp73 + tmp149 * tmp57 + tmp153 * tmp31 + tmp155 );
      DFinvDFx( 1, 3 ) = tmp127 * tmp160 + tmp87 * ( tmp148 * tmp77 + tmp161 + tmp31 * tmp91 + tmp65 * tmp68 );
      DFinvDFx( 1, 4 ) = tmp128 * tmp160 + tmp87 * ( tmp118 * tmp65 + tmp148 * tmp84 + tmp149 * tmp76 + tmp156 );
      DFinvDFx( 1, 5 ) = tmp131 * tmp160 + tmp87 * ( tmp118 * tmp30 + tmp141 * tmp84 + tmp150 * tmp78 + tmp69 * tmp80 );
      DFinvDFx( 2, 3 ) = tmp134 * tmp160 + tmp87 * ( tmp144 * tmp75 + tmp151 * tmp78 + tmp45 * tmp57 + tmp69 * tmp81 );
      DFinvDFx( 2, 4 ) = tmp136 * tmp160 + tmp87 * ( tmp118 * tmp47 + tmp140 * tmp84 + tmp153 * tmp76 + tmp73 * tmp91 );
      DFinvDFx( 2, 5 ) = tmp137 * tmp160 + tmp87 * ( tmp118 * tmp57 + tmp144 * tmp84 + tmp152 * tmp78 + tmp69 * tmp72 );
      DFinvDFx( 0, 6 ) = tmp171 * tmp88 + tmp87 * ( tmp142 * tmp65 + tmp147 * tmp31 + tmp164 + tmp165 * tmp57 );
      DFinvDFx( 0, 7 ) = tmp124 * tmp171 + tmp87 * ( tmp148 * tmp80 + tmp149 * tmp65 + tmp167 * tmp81 + tmp168 );
      DFinvDFx( 0, 8 ) = tmp125 * tmp171 + tmp87 * ( tmp142 * tmp79 + tmp150 * tmp31 + tmp163 * tmp73 + tmp166 * tmp57 );
      DFinvDFx( 1, 6 ) = tmp127 * tmp171 + tmp87 * ( tmp163 * tmp76 + tmp165 * tmp77 + tmp31 * tmp63 + tmp51 * tmp65 );
      DFinvDFx( 1, 7 ) = tmp128 * tmp171 + tmp87 * ( tmp113 * tmp65 + tmp165 * tmp84 + tmp166 * tmp76 + tmp63 * tmp79 );
      DFinvDFx( 1, 8 ) = tmp131 * tmp171 + tmp87 * ( tmp113 * tmp30 + tmp162 * tmp84 + tmp167 * tmp78 + tmp52 * tmp80 );
      DFinvDFx( 2, 6 ) = tmp134 * tmp171 + tmp87 * ( tmp142 * tmp75 + tmp148 * tmp78 + tmp52 * tmp81 + tmp57 * tmp62 );
      DFinvDFx( 2, 7 ) = tmp136 * tmp171 + tmp87 * ( tmp147 * tmp84 + tmp150 * tmp76 + tmp170 + tmp63 * tmp73 );
      DFinvDFx( 2, 8 ) = tmp137 * tmp171 + tmp87 * ( tmp113 * tmp57 + tmp142 * tmp84 + tmp149 * tmp78 + tmp169 );
   }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override
   {
      sendBuffer << Type::ICOSAHEDRAL_SHELL << rayVertex_ << refVertex_ << thrVertex_ << forVertex_ << radRefVertex_
                 << radRayVertex_ << prismNormal_;
   }

   static void setMap( SetupPrimitiveStorage& setupStorage )
   {
      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(), std::make_shared< IcosahedralShellMap >( cell, setupStorage ) );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(), std::make_shared< IcosahedralShellMap >( face, setupStorage ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(), std::make_shared< IcosahedralShellMap >( edge, setupStorage ) );
      }
   }

   const Point3D& rayVertex() const { return rayVertex_; }
   const Point3D& refVertex() const { return refVertex_; }
   const Point3D& thrVertex() const { return thrVertex_; }
   const Point3D& forVertex() const { return forVertex_; }

   real_t radRefVertex() const { return radRefVertex_; }
   real_t radRayVertex() const { return radRayVertex_; }

   /// Enumeration class for classifying tetrahedra
   enum class TetType
   {
      TET_INWARDS,
      TET_OUTWARDS,
      TET_SKEW
   };

   /// method for classifying the vertices of the macro tetrahedron
   static TetType classifyTet( const Cell& cell, std::array< real_t, 4 >& radius )
   {
      const std::array< Point3D, 4 >& coords = cell.getCoordinates();
      real_t                          innerRad, outerRad;
      innerRad = std::numeric_limits< real_t >::max();
      outerRad = 0.0;

      for ( uint_t k = 0; k < 4; k++ )
      {
         radius[k] = std::sqrt( coords[k].squaredNorm() );
         innerRad  = radius[k] < innerRad ? radius[k] : innerRad;
         outerRad  = radius[k] > outerRad ? radius[k] : outerRad;
      }

      SHELL_MAP_LOG( "outer radius = " << outerRad );
      SHELL_MAP_LOG( "inner radius = " << innerRad );

      uint_t nOuterNodes = 0;
      real_t eps         = ( outerRad - innerRad ) * real_c( 0.001 );
      for ( uint_t k = 0; k < 4; k++ )
      {
         SHELL_MAP_LOG( "radius[" << k << "] = " << radius[k] );
         nOuterNodes += ( outerRad - radius[k] ) < eps ? 1 : 0;
      }

      SHELL_MAP_LOG( "classifyTet: nOuterNodes = " << nOuterNodes );
      TetType thisTetType;
      switch ( nOuterNodes )
      {
      case 1:
         thisTetType = TetType::TET_OUTWARDS;
         SHELL_MAP_LOG( " -> TET_OUTWARDS" );
         break;
      case 2:
         thisTetType = TetType::TET_SKEW;
         SHELL_MAP_LOG( " -> TET_SKEW" );
         break;
      case 3:
         thisTetType = TetType::TET_INWARDS;
         SHELL_MAP_LOG( " -> TET_INWARDS" );
         break;
      default:
         WALBERLA_ABORT( "Houston we have a problem! Cannot classify macro tetrahedron!" );
      }

      return thisTetType;
   }

 private:
   /// \name Classified vertices of macro tetrahedron
   ///
   /// Each macro tetrahedron can be classified into one of three categories:
   /// - TET_OUTWARDS: one vertex on outer layer, three on inner layer
   /// - TET_INWARDS: one vertex on inner layer, three on outer layer
   /// - TET_SKEW: two vertices on inner layer, two on outer layer
   /// For the first two cases the single vertex and one of the other three are lying on a ray coming from the
   /// origin. The single one will become the refVertex_, the other one the rayVertex_. The remaining two
   /// two vertices become thrVertex_ and forVertex_.
   /// In the case of TET_SKEW a neighbouring tetrahedron from the same prism that is non-skew is found and
   /// its vertices used instead.
   ///@{
   Point3D refVertex_;
   Point3D rayVertex_;
   Point3D thrVertex_;
   Point3D forVertex_;
   ///@}

   /// tolerance for comparing numerical values for equality
   const real_t tol = real_c( std::is_same< real_t, double >() ? 1e-14 : 1e-7 );

   /// distance from origin of vertex rayVertex_
   real_t radRefVertex_;

   /// distance from origin of vertex refVertex_
   real_t radRayVertex_;

   /// normal of prism planes
   Point3D prismNormal_;

   void classifyVertices( const Cell& cell, const SetupPrimitiveStorage& storage )
   {
      const std::array< Point3D, 4 >& coords = cell.getCoordinates();

      real_t aux;
      bool   pairFound = false;

      uint_t idxRefVertex = 0, idxRayVertex = 0, idxThrVertex = 0, idxForVertex = 0;

      SHELL_MAP_LOG( "macro-vertex 0 = " << coords[0] );
      SHELL_MAP_LOG( "macro-vertex 1 = " << coords[1] );
      SHELL_MAP_LOG( "macro-vertex 2 = " << coords[2] );
      SHELL_MAP_LOG( "macro-vertex 3 = " << coords[3] );

      // determine type of macro-tet
      std::array< real_t, 4 > radius;
      TetType                 thisTetType = classifyTet( cell, radius );

      // determine the two vertices lying on a radial ray
      for ( uint_t k = 0; k < 4 && !pairFound; k++ )
      {
         for ( uint_t j = k + 1; j < 4 && !pairFound; j++ )
         {
            SHELL_MAP_LOG( "Testing cross product v[" << k << "] x v[" << j << "]" );

            // x-component of cross-product
            aux = coords[k][1] * coords[j][2] - coords[k][2] * coords[j][1];
            if ( std::abs( aux ) / ( radius[k] * radius[j] ) > tol )
               continue;

            // y-component of cross-product
            aux = coords[k][2] * coords[j][0] - coords[k][0] * coords[j][2];
            if ( std::abs( aux ) / ( radius[k] * radius[j] ) > tol )
               continue;

            // z-component of cross-product
            aux = coords[k][0] * coords[j][1] - coords[k][1] * coords[j][0];
            if ( std::abs( aux ) / ( radius[k] * radius[j] ) > tol )
               continue;

            // still here, so we found the vertex pair (ordering does not matter)
            idxRefVertex = k;
            idxRayVertex = j;
            pairFound    = true;
         }
      }

      if ( !pairFound )
      {
         WALBERLA_ABORT( "Error in finding vertex pair on radial ray!!!" );
      }

      // ------------------------------------------------------------------
      //  for skew tets we have a problem, so we need to find another tet
      //  from the same prims that is non-skew
      // ------------------------------------------------------------------
      if ( thisTetType == TetType::TET_SKEW )
      {
         PrimitiveID altCellID = findNonSkewTetInPrism( cell, storage, idxRefVertex, idxRayVertex );
         classifyVertices( *storage.getCell( altCellID ), storage );
      }

      else
      {
         // remember the indices of the two remaining nodes
         for ( uint_t k = 0; k < 4; k++ )
         {
            if ( k != idxRefVertex && k != idxRayVertex )
            {
               idxThrVertex = k;
            }
         }
         for ( uint_t k = 0; k < 4; k++ )
         {
            if ( k != idxRefVertex && k != idxRayVertex && k != idxThrVertex )
            {
               idxForVertex = k;
            }
         }

         // now sort ref and ray vertices depending on tet type
         switch ( thisTetType )
         {
         case TetType::TET_OUTWARDS:
            if ( radius[idxRefVertex] < radius[idxRayVertex] )
            {
               uint_t swp   = idxRayVertex;
               idxRayVertex = idxRefVertex;
               idxRefVertex = swp;
            }
            break;

         case TetType::TET_INWARDS:
            if ( radius[idxRefVertex] > radius[idxRayVertex] )
            {
               uint_t swp   = idxRayVertex;
               idxRayVertex = idxRefVertex;
               idxRefVertex = swp;
            }
            break;

         default:
            WALBERLA_ABORT( "We should not have ended up in this case!" );
         }

         // store vertex coordinates and radii
         refVertex_ = coords[idxRefVertex];
         rayVertex_ = coords[idxRayVertex];
         thrVertex_ = coords[idxThrVertex];
         forVertex_ = coords[idxForVertex];

         radRefVertex_ = std::sqrt( refVertex_.squaredNorm() );
         radRayVertex_ = std::sqrt( rayVertex_.squaredNorm() );

         SHELL_MAP_LOG( "refVertex = " << refVertex_ );
         SHELL_MAP_LOG( "rayVertex = " << rayVertex_ );
         SHELL_MAP_LOG( "thrVertex = " << thrVertex_ );
         SHELL_MAP_LOG( "forVertex = " << forVertex_ );

         // calculate normal of prism
         // ray, thr, and for are the three vertices spanning either the outer or inner plane
         prismNormal_ = crossProduct( thrVertex_ - rayVertex_, forVertex_ - rayVertex_ );
         prismNormal_ /= prismNormal_.norm();
      }
   }

   /// For a "skew" tetrahedron find a non-skew neighbour from the same prism
   PrimitiveID
       findNonSkewTetInPrism( const Cell& cell, const SetupPrimitiveStorage& storage, uint_t& idxRefVertex, uint_t& idxRayVertex )
   {
      std::vector< PrimitiveID > verts;
      cell.getNeighborVertices( verts );

      PrimitiveID rayNode = verts[idxRayVertex];
      PrimitiveID refNode = verts[idxRefVertex];

      WALBERLA_ASSERT_EQUAL( idxRayVertex, cell.getLocalVertexID( rayNode ) );
      WALBERLA_ASSERT_EQUAL( idxRefVertex, cell.getLocalVertexID( refNode ) );

      std::vector< PrimitiveID > nbrFaces;
      cell.getNeighborFaces( nbrFaces );
      SHELL_MAP_LOG( "Skew cell has " << nbrFaces.size() << " face neighbours" );

      for ( uint_t k = 0; k < nbrFaces.size(); k++ )
      {
         const std::shared_ptr< const Face > candidate = storage.getFace( nbrFaces[k] );
         candidate->getNeighborVertices( verts );

         // check that not both vertices are part of the face
         auto itRef = std::find( verts.begin(), verts.end(), refNode );
         auto itRay = std::find( verts.begin(), verts.end(), rayNode );
         if ( !( itRef != verts.end() && itRay != verts.end() ) )
         {
            SHELL_MAP_LOG( "Found fitting face: id = " << nbrFaces[k] );

            // now select the correct cell
            std::vector< PrimitiveID > nbrCells;
            candidate->getNeighborCells( nbrCells );
            SHELL_MAP_LOG( "  Candidate face has " << nbrCells.size() << " cell neighbours" );
            for ( uint_t j = 0; j < nbrCells.size(); j++ )
            {
               if ( nbrCells[j] != cell.getID() )
               {
                  SHELL_MAP_LOG( "    Cell we need has ID = " << nbrCells[j] );
                  return nbrCells[j];
               }
            }
         }
      }
      WALBERLA_ABORT( "Cound not findNonSkewTetInPrism()" );
   }
};

} // namespace hyteg
