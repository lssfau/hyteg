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

#define ANNULUS_MAP_LOG( MSG )
// #define ANNULUS_MAP_LOG( MSG ) WALBERLA_LOG_INFO_ON_ROOT( MSG );

namespace hyteg {

/// Class providing geometry mapping for an annulus
///
/// This class takes an annulus mesh generated using MeshInfo::meshAnnulus()
/// with the CRISS or CROSS flavour and maps nodes on refined mesh levels to
/// corresponding radial layers. The mapping is performed with local
/// information from the individual macro face (as we did in 3D in HHG).
/// Note that the resulting mesh is different from the one we get using
/// MeshInfo::meshRectangle() in combination with the PolarCoordsMap.
class AnnulusAlignedMap : public GeometryMap
{
 public:
   AnnulusAlignedMap( const Face& face )
   {
      ANNULUS_MAP_LOG( "---------------------------------------------" );
      ANNULUS_MAP_LOG( "Initialising Annulus map for faceID: " << face.getID() );
      ANNULUS_MAP_LOG( "---------------------------------------------" );

      classifyVertices( face.getCoordinates() );
   }

   AnnulusAlignedMap( const Edge& edge, const SetupPrimitiveStorage& storage )
   {
      ANNULUS_MAP_LOG( "---------------------------------------------" );
      ANNULUS_MAP_LOG( "Initialising Annulus map for edgeID: " << edge.getID() );
      ANNULUS_MAP_LOG( "---------------------------------------------" );

      std::vector< PrimitiveID > neighborFaces;
      edge.getNeighborFaces( neighborFaces );
      WALBERLA_ASSERT_GREATER( neighborFaces.size(), 0 );
      const Face& face = *storage.getFace( neighborFaces[0] );
      classifyVertices( face.getCoordinates() );
   }

   AnnulusAlignedMap( walberla::mpi::RecvBuffer& recvBuffer )
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

      recvBuffer >> radRefVertex_;
      recvBuffer >> radRayVertex_;
   }

   void evalF( const Point3D& xold, Point3D& xnew ) const override
   {
      real_t P_x = xold[0];
      real_t P_y = xold[1];
      real_t A_x = A[0];
      real_t A_y = A[1];
      real_t B_x = B[0];
      real_t B_y = B[1];
      real_t C_x = C[0];
      real_t C_y = C[1];
      real_t D_x = D[0];
      real_t D_y = D[1];

      real_t MtM_inv_Mt_0 = MtM_inv_Mt( 0 );
      real_t MtM_inv_Mt_1 = MtM_inv_Mt( 1 );

      real_t tmp0  = A_x - D_x;
      real_t tmp1  = -tmp0;
      real_t tmp2  = A_x - B_x;
      real_t tmp3  = -tmp2;
      real_t tmp4  = A_y - C_y;
      real_t tmp5  = A_x - C_x;
      real_t tmp6  = -tmp5;
      real_t tmp7  = A_y - B_y;
      real_t tmp8  = -tmp7;
      real_t tmp9  = A_x * tmp4 - P_x * tmp4;
      real_t tmp10 = ( A_y * tmp6 - P_y * tmp6 + tmp9 ) / ( tmp3 * tmp4 + tmp6 * tmp8 );
      real_t tmp11 = MtM_inv_Mt_0 * ( -A_x + P_x + tmp10 * tmp3 ) + MtM_inv_Mt_1 * ( -A_y + P_y + tmp10 * tmp8 );
      real_t tmp12 = A_x + tmp1 * tmp11;
      real_t tmp13 = ( -A_y * tmp5 + P_y * tmp5 + tmp9 ) / ( -tmp2 * tmp4 + tmp5 * tmp7 );
      real_t tmp14 = A_x - P_x + tmp13 * tmp2;
      real_t tmp15 = A_y - P_y + tmp13 * tmp7;
      real_t tmp16 = -MtM_inv_Mt_0 * tmp14 - MtM_inv_Mt_1 * tmp15;
      real_t tmp17 = A_y - D_y;
      real_t tmp18 = -MtM_inv_Mt_0 * tmp14 - MtM_inv_Mt_1 * tmp15;
      real_t tmp19 = P_x * tmp17;
      real_t tmp20 = ( -P_y * tmp0 + tmp19 ) / ( tmp0 * ( A_y - tmp17 * tmp18 ) - tmp17 * ( A_x - tmp0 * tmp18 ) );
      real_t tmp21 = ( -A_x * tmp17 - A_y * tmp1 + P_x * tmp17 + P_y * tmp1 ) / ( tmp1 * tmp8 + tmp17 * tmp3 );
      real_t tmp22 = A_y - tmp11 * tmp17;
      real_t tmp23 =
          ( P_y * tmp1 + tmp19 ) * sqrt( pow( fabs( A_x + tmp21 * tmp3 ), 2 ) + pow( fabs( A_y + tmp21 * tmp8 ), 2 ) ) /
          ( ( tmp1 * tmp22 + tmp12 * tmp17 ) *
            sqrt( pow( fabs( tmp20 * ( A_x - tmp0 * tmp16 ) ), 2 ) + pow( fabs( tmp20 * ( A_y - tmp16 * tmp17 ) ), 2 ) ) );
      xnew[0] = tmp12 * tmp23;
      xnew[1] = tmp22 * tmp23;

      ANNULUS_MAP_LOG( "Mapped: " << xold << " --> " << xnew );
   }

   void evalFinv( const Point3D& xPhys, Point3D& xComp ) const override
   {
      real_t P_tilde_x = xPhys[0];
      real_t P_tilde_y = xPhys[1];
      real_t A_x       = A[0];
      real_t A_y       = A[1];
      real_t B_x       = B[0];
      real_t B_y       = B[1];
      real_t C_x       = C[0];
      real_t C_y       = C[1];
      real_t D_x       = D[0];
      real_t D_y       = D[1];

      real_t NtN_inv_Nt_0 = NtN_inv_Nt( 0 );
      real_t NtN_inv_Nt_1 = NtN_inv_Nt( 1 );

      real_t tmp0  = -A_x + B_x;
      real_t tmp1  = A_y - D_y;
      real_t tmp2  = -A_x + D_x;
      real_t tmp3  = -A_y + B_y;
      real_t tmp4  = A_x * tmp1;
      real_t tmp5  = one_over_A_norm * sqrt( pow( P_tilde_x, 2 ) + pow( P_tilde_y, 2 ) );
      real_t tmp6  = A_y * tmp2;
      real_t tmp7  = tmp4 + tmp6;
      real_t tmp8  = 1.0 / ( P_tilde_x * tmp1 + P_tilde_y * tmp2 );
      real_t tmp9  = NtN_inv_Nt_0 * ( -A_x + P_tilde_x * tmp7 * tmp8 ) + NtN_inv_Nt_1 * ( -A_y + P_tilde_y * tmp7 * tmp8 );
      real_t tmp10 = A_y + tmp9 * ( -A_y + C_y );
      real_t tmp11 = A_x + tmp9 * ( -A_x + C_x );
      real_t tmp12 = ( -tmp1 * tmp11 - tmp10 * tmp2 + tmp4 * tmp5 + tmp5 * tmp6 ) / ( tmp0 * tmp1 + tmp2 * tmp3 );
      xComp[0]     = tmp0 * tmp12 + tmp11;
      xComp[1]     = tmp10 + tmp12 * tmp3;
   }

   void evalDF( const Point3D& x, Matrix2r& DFx ) const override
   {
      Point3D xnew;
      evalF( x, xnew );
      if ( ( x - xnew ).norm() < 1e-14 )
      {
         DFx.setIdentity();
         return;
      }

      real_t P_x = x[0];
      real_t P_y = x[1];
      real_t A_x = A[0];
      real_t A_y = A[1];
      real_t B_x = B[0];
      real_t B_y = B[1];
      real_t C_x = C[0];
      real_t C_y = C[1];
      real_t D_x = D[0];
      real_t D_y = D[1];

      real_t MtM_inv_Mt_0 = MtM_inv_Mt( 0 );
      real_t MtM_inv_Mt_1 = MtM_inv_Mt( 1 );

      real_t tmp0  = A_x - D_x;
      real_t tmp1  = -tmp0;
      real_t tmp2  = A_y - B_y;
      real_t tmp3  = -tmp2;
      real_t tmp4  = A_x - B_x;
      real_t tmp5  = -tmp4;
      real_t tmp6  = A_y - C_y;
      real_t tmp7  = A_x - C_x;
      real_t tmp8  = -tmp7;
      real_t tmp9  = 1.0 / ( tmp3 * tmp8 + tmp5 * tmp6 );
      real_t tmp10 = -tmp6;
      real_t tmp11 = tmp10 * tmp9;
      real_t tmp12 = MtM_inv_Mt_0 * ( tmp11 * tmp5 + 1 ) + MtM_inv_Mt_1 * tmp11 * tmp3;
      real_t tmp13 = tmp1 * tmp12;
      real_t tmp14 = A_y - D_y;
      real_t tmp15 = P_x * tmp14;
      real_t tmp16 = P_y * tmp1 + tmp15;
      real_t tmp17 = -tmp2 * tmp7 + tmp4 * tmp6;
      real_t tmp18 = -1 / tmp17;
      real_t tmp19 = A_x * tmp6 - P_x * tmp6;
      real_t tmp20 = -A_y * tmp7 + P_y * tmp7 + tmp19;
      real_t tmp21 = tmp18 * tmp20;
      real_t tmp22 = A_x - P_x + tmp21 * tmp4;
      real_t tmp23 = A_y - P_y + tmp2 * tmp21;
      real_t tmp24 = MtM_inv_Mt_0 * tmp22 + MtM_inv_Mt_1 * tmp23;
      real_t tmp25 = -tmp24;
      real_t tmp26 = A_x - tmp0 * tmp25;
      real_t tmp27 = -MtM_inv_Mt_0 * tmp22 - MtM_inv_Mt_1 * tmp23;
      real_t tmp28 = 1.0 / ( tmp0 * ( A_y - tmp14 * tmp27 ) - tmp14 * ( A_x - tmp0 * tmp27 ) );
      real_t tmp29 = -P_y * tmp0 + tmp15;
      real_t tmp30 = tmp28 * tmp29;
      real_t tmp31 = fabs( tmp26 * tmp30 );
      real_t tmp32 = A_y - tmp14 * tmp25;
      real_t tmp33 = fabs( tmp30 * tmp32 );
      real_t tmp34 = pow( tmp31, 2 ) + pow( tmp33, 2 );
      real_t tmp35 = pow( tmp34, -1.0 / 2.0 );
      real_t tmp36 = tmp14 * tmp5;
      real_t tmp37 = tmp1 * tmp3;
      real_t tmp38 = 1.0 / ( tmp36 + tmp37 );
      real_t tmp39 = tmp38 * ( -A_x * tmp14 - A_y * tmp1 + P_x * tmp14 + P_y * tmp1 );
      real_t tmp40 = A_x + tmp39 * tmp5;
      real_t tmp41 = fabs( tmp40 );
      real_t tmp42 = A_y + tmp3 * tmp39;
      real_t tmp43 = fabs( tmp42 );
      real_t tmp44 = sqrt( pow( tmp41, 2 ) + pow( tmp43, 2 ) );
      real_t tmp45 = -tmp14;
      real_t tmp46 = tmp9 * ( A_y * tmp8 - P_y * tmp8 + tmp19 );
      real_t tmp47 = -A_x + P_x;
      real_t tmp48 = -A_y + P_y;
      real_t tmp49 = MtM_inv_Mt_0 * ( tmp46 * tmp5 + tmp47 ) + MtM_inv_Mt_1 * ( tmp3 * tmp46 + tmp48 );
      real_t tmp50 = A_y + tmp45 * tmp49;
      real_t tmp51 = tmp1 * tmp50;
      real_t tmp52 = A_x + tmp1 * tmp49;
      real_t tmp53 = tmp14 * tmp52;
      real_t tmp54 = tmp51 + tmp53;
      real_t tmp55 = 1.0 / tmp54;
      real_t tmp56 = tmp44 * tmp55;
      real_t tmp57 = tmp35 * tmp56;
      real_t tmp58 = tmp16 * tmp57;
      real_t tmp59 = -tmp13 * tmp14 - tmp13 * tmp45;
      real_t tmp60 = tmp16 * tmp52;
      real_t tmp61 = tmp35 * tmp60;
      real_t tmp62 = tmp44 / pow( tmp54, 2 );
      real_t tmp63 = tmp61 * tmp62;
      real_t tmp64 = tmp38 * tmp41 * ( ( ( tmp40 ) > 0 ) - ( ( tmp40 ) < 0 ) );
      real_t tmp65 = tmp38 * tmp43 * ( ( ( tmp42 ) > 0 ) - ( ( tmp42 ) < 0 ) );
      real_t tmp66 = tmp14 * tmp3 * tmp65 + tmp36 * tmp64;
      real_t tmp67 = tmp55 / tmp44;
      real_t tmp68 = tmp61 * tmp67;
      real_t tmp69 = tmp10 * tmp18;
      real_t tmp70 = tmp30 * ( -MtM_inv_Mt_0 * ( tmp4 * tmp69 - 1 ) - MtM_inv_Mt_1 * tmp2 * tmp69 );
      real_t tmp71 = tmp20 / tmp17;
      real_t tmp72 = MtM_inv_Mt_0 * ( -tmp4 * tmp71 - tmp47 ) + MtM_inv_Mt_1 * ( -tmp2 * tmp71 - tmp48 );
      real_t tmp73 = A_x + tmp0 * tmp72;
      real_t tmp74 = tmp0 * ( A_y + tmp14 * tmp24 ) - tmp14 * ( A_x + tmp0 * tmp24 );
      real_t tmp75 = tmp29 / tmp74;
      real_t tmp76 = tmp28 * tmp74;
      real_t tmp77 = tmp26 * tmp31 * tmp76 * ( ( ( tmp73 * tmp75 ) > 0 ) - ( ( tmp73 * tmp75 ) < 0 ) ) / tmp73;
      real_t tmp78 = A_y + tmp14 * tmp72;
      real_t tmp79 = tmp32 * tmp33 * tmp76 * ( ( ( tmp75 * tmp78 ) > 0 ) - ( ( tmp75 * tmp78 ) < 0 ) ) / tmp78;
      real_t tmp80 = -tmp77 * ( -tmp0 * tmp70 + tmp14 * tmp26 * tmp28 ) - tmp79 * ( tmp14 * tmp28 * tmp32 - tmp14 * tmp70 );
      real_t tmp81 = tmp56 / pow( tmp34, 3.0 / 2.0 );
      real_t tmp82 = tmp60 * tmp81;
      real_t tmp83 = tmp7 * tmp9;
      real_t tmp84 = MtM_inv_Mt_0 * tmp5 * tmp83 + MtM_inv_Mt_1 * ( tmp3 * tmp83 + 1 );
      real_t tmp85 = tmp1 * tmp84;
      real_t tmp86 = -tmp14 * tmp85 - tmp45 * tmp85;
      real_t tmp87 = tmp1 * tmp5 * tmp64 + tmp37 * tmp65;
      real_t tmp88 = tmp30 * ( -MtM_inv_Mt_0 * tmp18 * tmp4 * tmp7 - MtM_inv_Mt_1 * ( tmp18 * tmp2 * tmp7 - 1 ) );
      real_t tmp89 = tmp1 * tmp28;
      real_t tmp90 = -tmp77 * ( -tmp0 * tmp88 + tmp26 * tmp89 ) - tmp79 * ( -tmp14 * tmp88 + tmp32 * tmp89 );
      real_t tmp91 = tmp45 * tmp58;
      real_t tmp92 = tmp16 * tmp50;
      real_t tmp93 = tmp35 * tmp92;
      real_t tmp94 = tmp62 * tmp93;
      real_t tmp95 = tmp67 * tmp93;
      real_t tmp96 = tmp81 * tmp92;
      DFx( 0, 0 )  = tmp13 * tmp58 + tmp53 * tmp57 + tmp59 * tmp63 + tmp66 * tmp68 + tmp80 * tmp82;
      DFx( 0, 1 )  = tmp1 * tmp52 * tmp57 + tmp58 * tmp85 + tmp63 * tmp86 + tmp68 * tmp87 + tmp82 * tmp90;
      DFx( 1, 0 )  = tmp12 * tmp91 + tmp14 * tmp50 * tmp57 + tmp59 * tmp94 + tmp66 * tmp95 + tmp80 * tmp96;
      DFx( 1, 1 )  = tmp51 * tmp57 + tmp84 * tmp91 + tmp86 * tmp94 + tmp87 * tmp95 + tmp90 * tmp96;
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const override
   {
      Matrix2r tmp;
      evalDF( x, tmp );
      DFinvx = tmp.inverse();
   }

   void evalDFinvDF( const Point3D& x, Matrixr< 2, 4 >& DFinvDFx ) const override final
   {
      WALBERLA_UNUSED( DFinvDFx );
      WALBERLA_ABORT( "Not implemented" );
   }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override
   {
      sendBuffer << Type::ANNULUS << rayVertex_ << refVertex_ << thrVertex_ << radRefVertex_ << radRayVertex_;
   }

   static void setMap( SetupPrimitiveStorage& setupStorage )
   {
      WALBERLA_LOG_WARNING_ON_ROOT(
          "The AnnulusAlignedMap (blending/geometry map) has not been well-tested for use with blending-capable operators.\n"
          "Remove this warning when this has been tested.\n"
          "Otherwise, it is safer to only use the map's evaluation function to initialize a parametric map (aka MicroMesh)." )

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(), std::make_shared< AnnulusAlignedMap >( face ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(), std::make_shared< AnnulusAlignedMap >( edge, setupStorage ) );
      }
   }

   const Point3D& rayVertex() const { return rayVertex_; }
   const Point3D& refVertex() const { return refVertex_; }
   const Point3D& thrVertex() const { return thrVertex_; }

   real_t radRefVertex() const { return radRefVertex_; }
   real_t radRayVertex() const { return radRayVertex_; }

 private:
   /// \name Classified vertices of macro triangle
   ///
   /// Each macro triangle of the annulus has two vertices which lie on a ray coming from the origin and
   /// two with the same distance from the origin. The vertex opposite to the edge formed by the latter
   /// two is stored as rayVertex_, the one on the ray with the rayVertex_ is stored as refVertex_.
   /// The third vertex then is stored as thrVertex_.
   ///@{
   Point3D rayVertex_;
   Point3D refVertex_;
   Point3D thrVertex_;
   ///@}

   /// distance from origin of vertex rayVertex_
   real_t radRefVertex_;

   /// distance from origin of vertex refVertex_
   real_t radRayVertex_;

   Point3D A;
   Point3D B;
   Point3D C;
   Point3D D;

   real_t one_over_A_norm;

   /// M = C - A
   /// This variable stores (M^T M)^{-1} M^T
   MatrixXr MtM_inv_Mt;

   /// N = D - A
   /// This variable stores (N^T N)^{-1} N^T
   MatrixXr NtN_inv_Nt;

   /// method for classifying the vertices of the macro triangle
   void classifyVertices( const std::array< Point3D, 3 >& coords )
   {
      std::array< real_t, 3 > radius{};
      for ( uint_t k = 0; k < 3; k++ )
      {
         radius[k] = std::sqrt( coords[k].squaredNorm() );
         ANNULUS_MAP_LOG( "Vertex " << k << ": (" // << std::scientific
                                    << coords[k][0] << ", " << coords[k][1] << ", " << coords[k][2] << ")\n"
                                    << " radius = " << radius[k] );
      }
      real_t cross01 = coords[0][0] * coords[1][1] - coords[0][1] * coords[1][0];
      real_t cross02 = coords[0][0] * coords[2][1] - coords[0][1] * coords[2][0];
      real_t cross12 = coords[1][0] * coords[2][1] - coords[1][1] * coords[2][0];

      ANNULUS_MAP_LOG( "r0 x r1 = " << std::showpos << std::scientific << cross01 );
      ANNULUS_MAP_LOG( "r0 x r2 = " << std::showpos << std::scientific << cross02 );
      ANNULUS_MAP_LOG( "r1 x r2 = " << std::showpos << std::scientific << cross12 );

      auto   dp     = std::is_same< real_t, double >();
      real_t tol    = real_c( dp ? 1e-14 : 1e-5 );
      uint_t intRay = 99;
      uint_t intRef = 99;

      // classify assuming we have a triangle pointing outwards from the origin
      if ( std::abs( cross01 ) < tol )
      {
         thrVertex_ = coords[2];
         if ( radius[0] < radius[1] )
         {
            intRay = 0;
            intRef = 1;
         }
         else
         {
            intRay = 1;
            intRef = 0;
         }
      }
      else if ( std::abs( cross02 ) < tol )
      {
         thrVertex_ = coords[1];
         if ( radius[0] < radius[2] )
         {
            intRay = 0;
            intRef = 2;
         }
         else
         {
            intRay = 2;
            intRef = 0;
         }
      }
      else if ( std::abs( cross12 ) < tol )
      {
         thrVertex_ = coords[0];
         if ( radius[1] < radius[2] )
         {
            intRay = 1;
            intRef = 2;
         }
         else
         {
            intRay = 2;
            intRef = 1;
         }
      }
      else
      {
         WALBERLA_ABORT( "Classification error in classifyVertices! Did you use CRISSCROSS maybe?" );
      }

      // swap classes in case we have a triangle pointing towards the origin
      ANNULUS_MAP_LOG( "Critical value = " << std::abs( coords[intRef].squaredNorm() - thrVertex_.squaredNorm() ) );
      if ( std::abs( coords[intRef].squaredNorm() - thrVertex_.squaredNorm() ) < tol )
      {
         ANNULUS_MAP_LOG( "Detected inward pointing triangle" );
         uint_t aux = intRay;
         intRay     = intRef;
         intRef     = aux;
      }
      else
      {
         ANNULUS_MAP_LOG( "Detected outward pointing triangle" );
      }

      rayVertex_ = coords[intRay];
      refVertex_ = coords[intRef];

      radRayVertex_ = radius[intRay];
      radRefVertex_ = radius[intRef];

      ANNULUS_MAP_LOG( "refVertex: (" << refVertex_[0] << ", " << refVertex_[1] << ")" );
      ANNULUS_MAP_LOG( "rayVertex: (" << rayVertex_[0] << ", " << rayVertex_[1] << ")" );
      ANNULUS_MAP_LOG( "thrVertex: (" << thrVertex_[0] << ", " << thrVertex_[1] << ")" );

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

      // C is the other vertex
      C = thrVertex_;

      // D is the other vertex on the outer boundary that spans the trapezoid with A
      // may be the same as C depending on the triangle type
      D = ( C / C.norm() ) * A.norm();

      one_over_A_norm = 1.0 / A.norm();

      MatrixXr M = C - A;
      MtM_inv_Mt = ( M.transpose() * M ).inverse() * M.transpose();

      MatrixXr N = D - A;
      NtN_inv_Nt = ( N.transpose() * N ).inverse() * N.transpose();
   }
};

} // namespace hyteg
