/*
* Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr, Andreas Burkhart, Nils Kohl.
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

#include "GeometryMap.hpp"
#include "IcosahedralShellMap.hpp"

// #define SHELL_MAP_LOG( STR ) WALBERLA_LOG_INFO_ON_ROOT( STR );
#define SHELL_MAP_LOG( STR )

namespace hyteg {

using walberla::real_c;
using TetType = IcosahedralShellMap::TetType;

/// Class providing geometry mapping for a facetted isosahedral shell
///
/// This geometry map provides a blending operation for a base mesh generated
/// with the inline meshSphericalShell generator. Geometric nodes on refined
/// hierarchy levels are projected onto spherical layers.
///
/// Compared to the IcosahedralShellMap, this slightly more expensive map also aligns
/// all nodes on beams through the origin.
/// (For a visual comparison in 2D have a look at the documentation of the AnnulusMap
/// and the AnnulusAlignedMap.)
class IcosahedralShellAlignedMap : public GeometryMap
{
 public:
   IcosahedralShellAlignedMap( const Cell& cell, const SetupPrimitiveStorage& storage )
   {
      SHELL_MAP_LOG( "---------------------------------------------" );
      SHELL_MAP_LOG( "Initialising Shell map for cellID: " << cell.getID() );
      SHELL_MAP_LOG( "---------------------------------------------" );
      classifyVertices( cell, storage );
   }

   IcosahedralShellAlignedMap( const Face& face, const SetupPrimitiveStorage& storage )
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

   IcosahedralShellAlignedMap( const Edge& edge, const SetupPrimitiveStorage& storage )
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

   IcosahedralShellAlignedMap( walberla::mpi::RecvBuffer& recvBuffer )
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
      real_t P_x = xold[0];
      real_t P_y = xold[1];
      real_t P_z = xold[2];

      real_t A_x = A[0];
      real_t A_y = A[1];
      real_t A_z = A[2];

      real_t B_x = B[0];
      real_t B_y = B[1];
      real_t B_z = B[2];

      real_t C1_x = C1[0];
      real_t C1_y = C1[1];
      real_t C1_z = C1[2];

      real_t C2_x = C2[0];
      real_t C2_y = C2[1];
      real_t C2_z = C2[2];

      real_t D1_x = D1[0];
      real_t D1_y = D1[1];
      real_t D1_z = D1[2];

      real_t D2_x = D2[0];
      real_t D2_y = D2[1];
      real_t D2_z = D2[2];

      real_t MtM_inv_Mt_0_0 = MtM_inv_Mt( 0, 0 );
      real_t MtM_inv_Mt_0_1 = MtM_inv_Mt( 0, 1 );
      real_t MtM_inv_Mt_0_2 = MtM_inv_Mt( 0, 2 );

      real_t MtM_inv_Mt_1_0 = MtM_inv_Mt( 1, 0 );
      real_t MtM_inv_Mt_1_1 = MtM_inv_Mt( 1, 1 );
      real_t MtM_inv_Mt_1_2 = MtM_inv_Mt( 1, 2 );

      real_t tmp0  = A_x - D1_x;
      real_t tmp1  = -tmp0;
      real_t tmp2  = -A_x;
      real_t tmp3  = A_x - B_x;
      real_t tmp4  = -tmp3;
      real_t tmp5  = A_y - C1_y;
      real_t tmp6  = -tmp5;
      real_t tmp7  = A_z - C2_z;
      real_t tmp8  = -tmp7;
      real_t tmp9  = A_y - C2_y;
      real_t tmp10 = -tmp9;
      real_t tmp11 = A_z - C1_z;
      real_t tmp12 = -tmp11;
      real_t tmp13 = -tmp10 * tmp12 + tmp6 * tmp8;
      real_t tmp14 = A_y - B_y;
      real_t tmp15 = -tmp14;
      real_t tmp16 = A_x - C1_x;
      real_t tmp17 = -tmp16;
      real_t tmp18 = A_x - C2_x;
      real_t tmp19 = -tmp18;
      real_t tmp20 = -tmp12 * tmp19 + tmp17 * tmp8;
      real_t tmp21 = -tmp20;
      real_t tmp22 = A_z - B_z;
      real_t tmp23 = -tmp22;
      real_t tmp24 = tmp10 * tmp17 - tmp19 * tmp6;
      real_t tmp25 = ( A_x * tmp13 + A_y * tmp21 + A_z * tmp24 - P_x * tmp13 - P_y * tmp21 - P_z * tmp24 ) /
                     ( tmp13 * tmp4 + tmp15 * tmp21 + tmp23 * tmp24 );
      real_t tmp26 = P_x + tmp2 + tmp25 * tmp4;
      real_t tmp27 = -A_y;
      real_t tmp28 = P_y + tmp15 * tmp25 + tmp27;
      real_t tmp29 = -A_z;
      real_t tmp30 = P_z + tmp23 * tmp25 + tmp29;
      real_t tmp31 = MtM_inv_Mt_0_0 * tmp26 + MtM_inv_Mt_0_1 * tmp28 + MtM_inv_Mt_0_2 * tmp30;
      real_t tmp32 = A_x - D2_x;
      real_t tmp33 = -tmp32;
      real_t tmp34 = MtM_inv_Mt_1_0 * tmp26 + MtM_inv_Mt_1_1 * tmp28 + MtM_inv_Mt_1_2 * tmp30;
      real_t tmp35 = A_x + tmp1 * tmp31 + tmp33 * tmp34;
      real_t tmp36 = -tmp11 * tmp9 + tmp5 * tmp7;
      real_t tmp37 = tmp16 * tmp9 - tmp18 * tmp5;
      real_t tmp38 = ( A_x * tmp36 - A_y * tmp20 + A_z * tmp37 - P_x * tmp36 + P_y * tmp20 - P_z * tmp37 ) /
                     ( tmp14 * tmp20 - tmp22 * tmp37 - tmp3 * tmp36 );
      real_t tmp39 = -A_x + P_x - tmp3 * tmp38;
      real_t tmp40 = -A_y + P_y - tmp14 * tmp38;
      real_t tmp41 = -A_z + P_z - tmp22 * tmp38;
      real_t tmp42 = MtM_inv_Mt_0_0 * tmp39 + MtM_inv_Mt_0_1 * tmp40 + MtM_inv_Mt_0_2 * tmp41;
      real_t tmp43 = MtM_inv_Mt_1_0 * tmp39 + MtM_inv_Mt_1_1 * tmp40 + MtM_inv_Mt_1_2 * tmp41;
      real_t tmp44 = tmp0 * tmp42 + tmp2 + tmp32 * tmp43;
      real_t tmp45 = A_y - D2_y;
      real_t tmp46 = A_y - D1_y;
      real_t tmp47 = tmp0 * tmp45 - tmp32 * tmp46;
      real_t tmp48 = A_z - D1_z;
      real_t tmp49 = A_z - D2_z;
      real_t tmp50 = tmp29 + tmp42 * tmp48 + tmp43 * tmp49;
      real_t tmp51 = -tmp45 * tmp48 + tmp46 * tmp49;
      real_t tmp52 = tmp0 * tmp49 - tmp32 * tmp48;
      real_t tmp53 = tmp27 + tmp42 * tmp46 + tmp43 * tmp45;
      real_t tmp54 = ( P_x * tmp51 - P_y * tmp52 + P_z * tmp47 ) / ( -tmp44 * tmp51 - tmp47 * tmp50 + tmp52 * tmp53 );
      real_t tmp55 = -tmp46;
      real_t tmp56 = -tmp49;
      real_t tmp57 = -tmp45;
      real_t tmp58 = -tmp48;
      real_t tmp59 = tmp55 * tmp56 - tmp57 * tmp58;
      real_t tmp60 = -tmp1 * tmp56 + tmp33 * tmp58;
      real_t tmp61 = tmp1 * tmp57 - tmp33 * tmp55;
      real_t tmp62 = ( -A_x * tmp59 - A_y * tmp60 - A_z * tmp61 + P_x * tmp59 + P_y * tmp60 + P_z * tmp61 ) /
                     ( tmp15 * tmp60 + tmp23 * tmp61 + tmp4 * tmp59 );
      real_t tmp63 = A_z + tmp31 * tmp58 + tmp34 * tmp56;
      real_t tmp64 = A_y + tmp31 * tmp55 + tmp34 * tmp57;
      real_t tmp65 =
          ( P_x * tmp59 + P_y * tmp60 + P_z * tmp61 ) *
          sqrt( pow( fabs( A_x + tmp4 * tmp62 ), 2 ) + pow( fabs( A_y + tmp15 * tmp62 ), 2 ) +
                pow( fabs( A_z + tmp23 * tmp62 ), 2 ) ) /
          ( ( tmp35 * tmp59 + tmp60 * tmp64 + tmp61 * tmp63 ) *
            sqrt( pow( fabs( tmp44 * tmp54 ), 2 ) + pow( fabs( tmp50 * tmp54 ), 2 ) + pow( fabs( tmp53 * tmp54 ), 2 ) ) );
      xnew[0] = tmp35 * tmp65;
      xnew[1] = tmp64 * tmp65;
      xnew[2] = tmp63 * tmp65;

      // SHELL_MAP_LOG( "Mapped: " << xold << " --> " << xnew );
   }

   void evalFinv( const Point3D& xPhys, Point3D& xComp ) const override
   {
      real_t P_tilde_x = xPhys[0];
      real_t P_tilde_y = xPhys[1];
      real_t P_tilde_z = xPhys[2];

      real_t A_x = A[0];
      real_t A_y = A[1];
      real_t A_z = A[2];

      real_t B_x = B[0];
      real_t B_y = B[1];
      real_t B_z = B[2];

      real_t C1_x = C1[0];
      real_t C1_y = C1[1];
      real_t C1_z = C1[2];

      real_t C2_x = C2[0];
      real_t C2_y = C2[1];
      real_t C2_z = C2[2];

      real_t D1_x = D1[0];
      real_t D1_y = D1[1];
      real_t D1_z = D1[2];

      real_t D2_x = D2[0];
      real_t D2_y = D2[1];
      real_t D2_z = D2[2];

      real_t NtN_inv_Nt_0_0 = NtN_inv_Nt( 0, 0 );
      real_t NtN_inv_Nt_0_1 = NtN_inv_Nt( 0, 1 );
      real_t NtN_inv_Nt_0_2 = NtN_inv_Nt( 0, 2 );

      real_t NtN_inv_Nt_1_0 = NtN_inv_Nt( 1, 0 );
      real_t NtN_inv_Nt_1_1 = NtN_inv_Nt( 1, 1 );
      real_t NtN_inv_Nt_1_2 = NtN_inv_Nt( 1, 2 );

      real_t tmp0  = -A_x + B_x;
      real_t tmp1  = -A_y + D1_y;
      real_t tmp2  = -A_z + D2_z;
      real_t tmp3  = -A_y + D2_y;
      real_t tmp4  = -A_z + D1_z;
      real_t tmp5  = tmp1 * tmp2 - tmp3 * tmp4;
      real_t tmp6  = -A_y + B_y;
      real_t tmp7  = -A_x + D1_x;
      real_t tmp8  = -A_x + D2_x;
      real_t tmp9  = -tmp2 * tmp7 + tmp4 * tmp8;
      real_t tmp10 = -A_z + B_z;
      real_t tmp11 = -tmp1 * tmp8 + tmp3 * tmp7;
      real_t tmp12 = A_x * tmp5;
      real_t tmp13 = one_over_A_norm * sqrt( pow( P_tilde_x, 2 ) + pow( P_tilde_y, 2 ) + pow( P_tilde_z, 2 ) );
      real_t tmp14 = A_y * tmp9;
      real_t tmp15 = A_z * tmp11;
      real_t tmp16 = tmp12 + tmp14 + tmp15;
      real_t tmp17 = 1.0 / ( P_tilde_x * tmp5 + P_tilde_y * tmp9 + P_tilde_z * tmp11 );
      real_t tmp18 = -A_x + P_tilde_x * tmp16 * tmp17;
      real_t tmp19 = -A_y + P_tilde_y * tmp16 * tmp17;
      real_t tmp20 = -A_z + P_tilde_z * tmp16 * tmp17;
      real_t tmp21 = NtN_inv_Nt_0_0 * tmp18 + NtN_inv_Nt_0_1 * tmp19 + NtN_inv_Nt_0_2 * tmp20;
      real_t tmp22 = NtN_inv_Nt_1_0 * tmp18 + NtN_inv_Nt_1_1 * tmp19 + NtN_inv_Nt_1_2 * tmp20;
      real_t tmp23 = A_z + tmp21 * ( -A_z + C1_z ) + tmp22 * ( -A_z + C2_z );
      real_t tmp24 = A_y + tmp21 * ( -A_y + C1_y ) + tmp22 * ( -A_y + C2_y );
      real_t tmp25 = A_x + tmp21 * ( -A_x + C1_x ) + tmp22 * ( -A_x + C2_x );
      real_t tmp26 = ( -tmp11 * tmp23 + tmp12 * tmp13 + tmp13 * tmp14 + tmp13 * tmp15 - tmp24 * tmp9 - tmp25 * tmp5 ) /
                     ( tmp0 * tmp5 + tmp10 * tmp11 + tmp6 * tmp9 );
      xComp[0] = tmp0 * tmp26 + tmp25;
      xComp[1] = tmp24 + tmp26 * tmp6;
      xComp[2] = tmp10 * tmp26 + tmp23;
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
      real_t P_x = x[0];
      real_t P_y = x[1];
      real_t P_z = x[2];

      real_t A_x = A[0];
      real_t A_y = A[1];
      real_t A_z = A[2];

      real_t B_x = B[0];
      real_t B_y = B[1];
      real_t B_z = B[2];

      real_t C1_x = C1[0];
      real_t C1_y = C1[1];
      real_t C1_z = C1[2];

      real_t C2_x = C2[0];
      real_t C2_y = C2[1];
      real_t C2_z = C2[2];

      real_t D1_x = D1[0];
      real_t D1_y = D1[1];
      real_t D1_z = D1[2];

      real_t D2_x = D2[0];
      real_t D2_y = D2[1];
      real_t D2_z = D2[2];

      real_t MtM_inv_Mt_0_0 = MtM_inv_Mt( 0, 0 );
      real_t MtM_inv_Mt_0_1 = MtM_inv_Mt( 0, 1 );
      real_t MtM_inv_Mt_0_2 = MtM_inv_Mt( 0, 2 );

      real_t MtM_inv_Mt_1_0 = MtM_inv_Mt( 1, 0 );
      real_t MtM_inv_Mt_1_1 = MtM_inv_Mt( 1, 1 );
      real_t MtM_inv_Mt_1_2 = MtM_inv_Mt( 1, 2 );

      real_t tmp0   = -A_z + D1_z;
      real_t tmp1   = -A_y + B_y;
      real_t tmp2   = -A_x + B_x;
      real_t tmp3   = -A_y + C1_y;
      real_t tmp4   = -A_z + C2_z;
      real_t tmp5   = -A_y + C2_y;
      real_t tmp6   = -A_z + C1_z;
      real_t tmp7   = tmp3 * tmp4 - tmp5 * tmp6;
      real_t tmp8   = -A_x + C1_x;
      real_t tmp9   = -A_x + C2_x;
      real_t tmp10  = tmp4 * tmp8 - tmp6 * tmp9;
      real_t tmp11  = -tmp10;
      real_t tmp12  = -A_z + B_z;
      real_t tmp13  = -tmp3 * tmp9 + tmp5 * tmp8;
      real_t tmp14  = 1.0 / ( tmp1 * tmp11 + tmp12 * tmp13 + tmp2 * tmp7 );
      real_t tmp15  = -tmp14 * tmp7;
      real_t tmp16  = tmp1 * tmp15;
      real_t tmp17  = tmp12 * tmp15;
      real_t tmp18  = tmp15 * tmp2 + 1;
      real_t tmp19  = MtM_inv_Mt_0_0 * tmp18 + MtM_inv_Mt_0_1 * tmp16 + MtM_inv_Mt_0_2 * tmp17;
      real_t tmp20  = -A_z + D2_z;
      real_t tmp21  = MtM_inv_Mt_1_0 * tmp18 + MtM_inv_Mt_1_1 * tmp16 + MtM_inv_Mt_1_2 * tmp17;
      real_t tmp22  = tmp0 * tmp19 + tmp20 * tmp21;
      real_t tmp23  = -A_x + D1_x;
      real_t tmp24  = -A_y + D2_y;
      real_t tmp25  = -A_x + D2_x;
      real_t tmp26  = -A_y + D1_y;
      real_t tmp27  = tmp23 * tmp24 - tmp25 * tmp26;
      real_t tmp28  = A_x * tmp7 + A_y * tmp11 + A_z * tmp13 - P_x * tmp7;
      real_t tmp29  = -P_y * tmp11 + tmp28;
      real_t tmp30  = tmp14 * ( -P_z * tmp13 + tmp29 );
      real_t tmp31  = -A_x + P_x;
      real_t tmp32  = tmp2 * tmp30 + tmp31;
      real_t tmp33  = -A_y;
      real_t tmp34  = P_y + tmp33;
      real_t tmp35  = tmp1 * tmp30 + tmp34;
      real_t tmp36  = -A_z;
      real_t tmp37  = P_z + tmp12 * tmp30 + tmp36;
      real_t tmp38  = MtM_inv_Mt_0_0 * tmp32 + MtM_inv_Mt_0_1 * tmp35 + MtM_inv_Mt_0_2 * tmp37;
      real_t tmp39  = MtM_inv_Mt_1_0 * tmp32 + MtM_inv_Mt_1_1 * tmp35 + MtM_inv_Mt_1_2 * tmp37;
      real_t tmp40  = A_z + tmp0 * tmp38 + tmp20 * tmp39;
      real_t tmp41  = tmp27 * tmp40;
      real_t tmp42  = tmp0 * tmp25 - tmp20 * tmp23;
      real_t tmp43  = A_y + tmp24 * tmp39 + tmp26 * tmp38;
      real_t tmp44  = tmp42 * tmp43;
      real_t tmp45  = -tmp0 * tmp24 + tmp20 * tmp26;
      real_t tmp46  = A_x + tmp23 * tmp38 + tmp25 * tmp39;
      real_t tmp47  = tmp45 * tmp46;
      real_t tmp48  = tmp41 + tmp44 + tmp47;
      real_t tmp49  = 1.0 / tmp48;
      real_t tmp50  = P_x * tmp45;
      real_t tmp51  = P_y * tmp42 + tmp50;
      real_t tmp52  = P_z * tmp27 + tmp51;
      real_t tmp53  = tmp49 * tmp52;
      real_t tmp54  = tmp40 * tmp49;
      real_t tmp55  = tmp19 * tmp26 + tmp21 * tmp24;
      real_t tmp56  = tmp19 * tmp23 + tmp21 * tmp25;
      real_t tmp57  = pow( tmp48, -2 );
      real_t tmp58  = tmp52 * tmp57;
      real_t tmp59  = tmp58 * ( -tmp22 * tmp27 - tmp42 * tmp55 - tmp45 * tmp56 );
      real_t tmp60  = tmp22 * tmp53 + tmp40 * tmp59 + tmp45 * tmp54;
      real_t tmp61  = 1.0 / ( A_x * tmp45 + A_y * tmp42 + A_z * tmp27 );
      real_t tmp62  = pow( tmp40, 2 ) * pow( tmp52, 2 ) * tmp57;
      real_t tmp63  = tmp41 * tmp53;
      real_t tmp64  = tmp51 + tmp63;
      real_t tmp65  = tmp52 * tmp54;
      real_t tmp66  = -tmp13 * tmp65;
      real_t tmp67  = tmp14 * ( tmp29 + tmp66 );
      real_t tmp68  = tmp2 * tmp67 + tmp31;
      real_t tmp69  = tmp1 * tmp67 + tmp34;
      real_t tmp70  = tmp36 + tmp65;
      real_t tmp71  = tmp12 * tmp67 + tmp70;
      real_t tmp72  = MtM_inv_Mt_0_0 * tmp68 + MtM_inv_Mt_0_1 * tmp69 + MtM_inv_Mt_0_2 * tmp71;
      real_t tmp73  = MtM_inv_Mt_1_0 * tmp68 + MtM_inv_Mt_1_1 * tmp69 + MtM_inv_Mt_1_2 * tmp71;
      real_t tmp74  = A_y + tmp24 * tmp73 + tmp26 * tmp72;
      real_t tmp75  = tmp42 * tmp74;
      real_t tmp76  = tmp27 * ( A_z + tmp0 * tmp72 + tmp20 * tmp73 ) + tmp45 * ( A_x + tmp23 * tmp72 + tmp25 * tmp73 ) + tmp75;
      real_t tmp77  = pow( tmp64, 2 ) * pow( tmp74, 2 ) / pow( tmp76, 2 );
      real_t tmp78  = tmp64 / tmp76;
      real_t tmp79  = tmp75 * tmp78;
      real_t tmp80  = tmp63 + tmp79;
      real_t tmp81  = tmp50 + tmp80;
      real_t tmp82  = tmp74 * tmp78;
      real_t tmp83  = tmp14 * ( -tmp11 * tmp82 + tmp28 + tmp66 );
      real_t tmp84  = tmp2 * tmp83 + tmp31;
      real_t tmp85  = tmp12 * tmp83 + tmp70;
      real_t tmp86  = tmp1 * tmp83 + tmp33 + tmp82;
      real_t tmp87  = MtM_inv_Mt_0_0 * tmp84 + MtM_inv_Mt_0_1 * tmp86 + MtM_inv_Mt_0_2 * tmp85;
      real_t tmp88  = MtM_inv_Mt_1_0 * tmp84 + MtM_inv_Mt_1_1 * tmp86 + MtM_inv_Mt_1_2 * tmp85;
      real_t tmp89  = A_x + tmp23 * tmp87 + tmp25 * tmp88;
      real_t tmp90  = tmp45 * tmp89;
      real_t tmp91  = tmp27 * ( A_z + tmp0 * tmp87 + tmp20 * tmp88 ) + tmp42 * ( A_y + tmp24 * tmp88 + tmp26 * tmp87 ) + tmp90;
      real_t tmp92  = pow( tmp81, 2 ) * pow( tmp89, 2 ) / pow( tmp91, 2 );
      real_t tmp93  = tmp62 + tmp77 + tmp92;
      real_t tmp94  = tmp61 / sqrt( tmp93 );
      real_t tmp95  = tmp27 * tmp94;
      real_t tmp96  = tmp81 / tmp91;
      real_t tmp97  = tmp89 * tmp96;
      real_t tmp98  = tmp90 * tmp96;
      real_t tmp99  = tmp80 + tmp98;
      real_t tmp100 = tmp61 * tmp99 / pow( tmp93, 3.0 / 2.0 );
      real_t tmp101 = tmp100 * tmp97;
      real_t tmp102 = -tmp101 * tmp65;
      real_t tmp103 = tmp102 + tmp95 * tmp97;
      real_t tmp104 = tmp43 * tmp49;
      real_t tmp105 = tmp104 * tmp45 + tmp43 * tmp59 + tmp53 * tmp55;
      real_t tmp106 = tmp42 * tmp94;
      real_t tmp107 = -tmp101 * tmp82;
      real_t tmp108 = tmp106 * tmp97 + tmp107;
      real_t tmp109 = tmp46 * tmp59 + tmp47 * tmp49 + tmp53 * tmp56;
      real_t tmp110 = tmp94 * tmp99;
      real_t tmp111 = -tmp100 * tmp92 + tmp110 + tmp94 * tmp98;
      real_t tmp112 = tmp10 * tmp14;
      real_t tmp113 = tmp112 * tmp2;
      real_t tmp114 = tmp112 * tmp12;
      real_t tmp115 = tmp1 * tmp112 + 1;
      real_t tmp116 = MtM_inv_Mt_0_0 * tmp113 + MtM_inv_Mt_0_1 * tmp115 + MtM_inv_Mt_0_2 * tmp114;
      real_t tmp117 = MtM_inv_Mt_1_0 * tmp113 + MtM_inv_Mt_1_1 * tmp115 + MtM_inv_Mt_1_2 * tmp114;
      real_t tmp118 = tmp0 * tmp116 + tmp117 * tmp20;
      real_t tmp119 = tmp116 * tmp26 + tmp117 * tmp24;
      real_t tmp120 = tmp116 * tmp23 + tmp117 * tmp25;
      real_t tmp121 = tmp58 * ( -tmp118 * tmp27 - tmp119 * tmp42 - tmp120 * tmp45 );
      real_t tmp122 = tmp118 * tmp53 + tmp121 * tmp40 + tmp42 * tmp54;
      real_t tmp123 = tmp119 * tmp53 + tmp121 * tmp43 + tmp44 * tmp49;
      real_t tmp124 = tmp46 * tmp49;
      real_t tmp125 = tmp120 * tmp53 + tmp121 * tmp46 + tmp124 * tmp42;
      real_t tmp126 = -tmp13 * tmp14;
      real_t tmp127 = tmp126 * tmp2;
      real_t tmp128 = tmp1 * tmp126;
      real_t tmp129 = tmp12 * tmp126 + 1;
      real_t tmp130 = MtM_inv_Mt_0_0 * tmp127 + MtM_inv_Mt_0_1 * tmp128 + MtM_inv_Mt_0_2 * tmp129;
      real_t tmp131 = MtM_inv_Mt_1_0 * tmp127 + MtM_inv_Mt_1_1 * tmp128 + MtM_inv_Mt_1_2 * tmp129;
      real_t tmp132 = tmp0 * tmp130 + tmp131 * tmp20;
      real_t tmp133 = tmp130 * tmp26 + tmp131 * tmp24;
      real_t tmp134 = tmp130 * tmp23 + tmp131 * tmp25;
      real_t tmp135 = tmp58 * ( -tmp132 * tmp27 - tmp133 * tmp42 - tmp134 * tmp45 );
      real_t tmp136 = tmp132 * tmp53 + tmp135 * tmp40 + tmp41 * tmp49;
      real_t tmp137 = tmp104 * tmp27 + tmp133 * tmp53 + tmp135 * tmp43;
      real_t tmp138 = tmp124 * tmp27 + tmp134 * tmp53 + tmp135 * tmp46;
      real_t tmp139 = -tmp100 * tmp65 * tmp82;
      real_t tmp140 = tmp139 + tmp82 * tmp95;
      real_t tmp141 = tmp45 * tmp94;
      real_t tmp142 = tmp107 + tmp141 * tmp82;
      real_t tmp143 = -tmp100 * tmp77 + tmp110 + tmp79 * tmp94;
      real_t tmp144 = tmp106 * tmp65 + tmp139;
      real_t tmp145 = tmp102 + tmp141 * tmp65;
      real_t tmp146 = -tmp100 * tmp62 + tmp110 + tmp63 * tmp94;
      DFx( 0, 0 )   = tmp103 * tmp60 + tmp105 * tmp108 + tmp109 * tmp111;
      DFx( 0, 1 )   = tmp103 * tmp122 + tmp108 * tmp123 + tmp111 * tmp125;
      DFx( 0, 2 )   = tmp103 * tmp136 + tmp108 * tmp137 + tmp111 * tmp138;
      DFx( 1, 0 )   = tmp105 * tmp143 + tmp109 * tmp142 + tmp140 * tmp60;
      DFx( 1, 1 )   = tmp122 * tmp140 + tmp123 * tmp143 + tmp125 * tmp142;
      DFx( 1, 2 )   = tmp136 * tmp140 + tmp137 * tmp143 + tmp138 * tmp142;
      DFx( 2, 0 )   = tmp105 * tmp144 + tmp109 * tmp145 + tmp146 * tmp60;
      DFx( 2, 1 )   = tmp122 * tmp146 + tmp123 * tmp144 + tmp125 * tmp145;
      DFx( 2, 2 )   = tmp136 * tmp146 + tmp137 * tmp144 + tmp138 * tmp145;

      return DFx.determinant();
   };

   void evalDFinv( const Point3D& x, Matrix3r& DFx ) const final override
   {
      Matrix3r tmp;
      evalDF( x, tmp );
      DFx = tmp.inverse();
   }

   void evalDFinvDF( const Point3D& x, Matrixr< 3, 9 >& DFinvDFx ) const override final
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFinvDFx );
      WALBERLA_ABORT( "Not implemented" );
   }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override
   {
      sendBuffer << Type::ICOSAHEDRAL_SHELL << rayVertex_ << refVertex_ << thrVertex_ << forVertex_ << radRefVertex_
                 << radRayVertex_ << prismNormal_;
   }

   static void setMap( SetupPrimitiveStorage& setupStorage )
   {
      WALBERLA_LOG_WARNING_ON_ROOT(
          "The IcosahedralShellAlignedMap (blending/geometry map) has not been well-tested for use with blending-capable operators.\n"
          "Remove this warning when this has been tested.\n"
          "Otherwise, it is safer to only use the map's evaluation function to initialize a parametric map (aka MicroMesh)." )

      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(), std::make_shared< IcosahedralShellAlignedMap >( cell, setupStorage ) );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(), std::make_shared< IcosahedralShellAlignedMap >( face, setupStorage ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(), std::make_shared< IcosahedralShellAlignedMap >( edge, setupStorage ) );
      }
   }

   const Point3D& rayVertex() const { return rayVertex_; }
   const Point3D& refVertex() const { return refVertex_; }
   const Point3D& thrVertex() const { return thrVertex_; }
   const Point3D& forVertex() const { return forVertex_; }

   real_t radRefVertex() const { return radRefVertex_; }
   real_t radRayVertex() const { return radRayVertex_; }

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

   Point3D A;
   Point3D B;
   Point3D C1;
   Point3D C2;
   Point3D D1;
   Point3D D2;

   real_t one_over_A_norm;

   /// M = (C_1 - A, C_2 - A)
   /// This variable stores (M^T M)^{-1} M^T
   MatrixXr MtM_inv_Mt;

   /// N =  (D_1 - A, D_2 - A)
   /// This variable stores (N^T N)^{-1} N^T
   MatrixXr NtN_inv_Nt;

   void classifyVertices( const Cell& cell, const SetupPrimitiveStorage& storage )
   {
      WALBERLA_UNUSED( storage );

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
      TetType                 thisTetType = IcosahedralShellMap::classifyTet( cell, radius );

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

         case TetType::TET_SKEW:
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

         // A must be the outer vertex that lies on the ray
         // B must be the other (inner) vertex on the ray
         if ( radRayVertex_ > radRefVertex_ )
         {
            A = rayVertex_;
            B = refVertex_;
         }
         else
         {
            B = rayVertex_;
            A = refVertex_;
         }

         // C1 and C2 are the other vertices
         C1 = thrVertex_;
         C2 = forVertex_;

         // D1, D2 are the other vertices on the outer boundary that span the prism with A
         // may be the same as C1 C2 depending on the tet type
         D1 = ( C1 / C1.norm() ) * A.norm();
         D2 = ( C2 / C2.norm() ) * A.norm();

         one_over_A_norm = 1.0 / A.norm();

         MatrixXr M( 3, 2 );
         M.col( 0 ) = C1 - A;
         M.col( 1 ) = C2 - A;
         MtM_inv_Mt = ( M.transpose() * M ).inverse() * M.transpose();

         MatrixXr N( 3, 2 );
         N.col( 0 ) = D1 - A;
         N.col( 1 ) = D2 - A;
         NtN_inv_Nt = ( N.transpose() * N ).inverse() * N.transpose();
      }
   }
};

} // namespace hyteg
