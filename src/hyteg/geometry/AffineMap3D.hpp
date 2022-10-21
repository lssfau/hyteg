/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Marcus Mohr.
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

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

using walberla::real_c;

/// The class implements a generic affine mapping in 3D
///
/// The affine mapping is characterised by a matrix \f$M\f$ and a vector \f$v\f$
/// and defined as
/// \f[
/// x \mapsto M x + v
/// \f]
class AffineMap3D : public GeometryMap
{
 public:
   AffineMap3D( const Matrix3r& mat, const Point3D& vec )
   : mat_( mat )
   , vec_( vec )
   {
      // precompute Jacobian determinant
      real_t tmp1 = +mat_( 0, 0 ) * mat_( 1, 1 ) * mat_( 2, 2 );
      real_t tmp2 = -mat_( 0, 0 ) * mat_( 2, 1 ) * mat_( 1, 2 );
      real_t tmp3 = -mat_( 1, 0 ) * mat_( 0, 1 ) * mat_( 2, 2 );
      real_t tmp4 = +mat_( 1, 0 ) * mat_( 2, 1 ) * mat_( 0, 2 );
      real_t tmp5 = +mat_( 2, 0 ) * mat_( 0, 1 ) * mat_( 1, 2 );
      real_t tmp6 = -mat_( 2, 0 ) * mat_( 1, 1 ) * mat_( 0, 2 );
      jacDet_     = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6;
   }

   AffineMap3D( walberla::mpi::RecvBuffer& recvBuffer )
   {
      for ( uint_t i = 0; i < 3; i++ )
      {
         for ( uint_t j = 0; j < 3; j++ )
         {
            recvBuffer >> mat_( i, j );
         }
      }
      recvBuffer >> vec_[0];
      recvBuffer >> vec_[1];
      recvBuffer >> vec_[2];
   }

   void evalF( const Point3D& xold, Point3D& xnew ) const override final
   {
      xnew[0] = mat_( 0, 0 ) * xold[0] + mat_( 0, 1 ) * xold[1] + mat_( 0, 2 ) * xold[2] + vec_[0];
      xnew[1] = mat_( 1, 0 ) * xold[0] + mat_( 1, 1 ) * xold[1] + mat_( 1, 2 ) * xold[2] + vec_[1];
      xnew[2] = mat_( 2, 0 ) * xold[0] + mat_( 2, 1 ) * xold[1] + mat_( 2, 2 ) * xold[2] + vec_[2];
   }

   void evalFinv( const Point3D& xPhys, Point3D& xComp ) const
   {
      real_t tmp0 = -vec_[0] + xPhys[0];
      real_t tmp1 = -vec_[1] + xPhys[1];
      real_t tmp2 = -vec_[2] + xPhys[2];

      xComp[0] = tmp0 * ( mat_( 1, 1 ) * mat_( 2, 2 ) - mat_( 2, 1 ) * mat_( 1, 2 ) ) +
                 tmp1 * ( -mat_( 0, 1 ) * mat_( 2, 2 ) + mat_( 2, 1 ) * mat_( 0, 2 ) ) +
                 tmp2 * ( mat_( 0, 1 ) * mat_( 1, 2 ) - mat_( 1, 1 ) * mat_( 0, 2 ) );
      xComp[1] = tmp0 * ( -mat_( 1, 0 ) * mat_( 2, 2 ) + mat_( 2, 0 ) * mat_( 1, 2 ) ) +
                 tmp1 * ( mat_( 0, 0 ) * mat_( 2, 2 ) - mat_( 2, 0 ) * mat_( 0, 2 ) ) +
                 tmp2 * ( -mat_( 0, 0 ) * mat_( 1, 2 ) + mat_( 1, 0 ) * mat_( 0, 2 ) );
      xComp[2] = tmp0 * ( mat_( 1, 0 ) * mat_( 2, 1 ) - mat_( 2, 0 ) * mat_( 1, 1 ) ) +
                 tmp1 * ( -mat_( 0, 0 ) * mat_( 2, 1 ) + mat_( 2, 0 ) * mat_( 0, 1 ) ) +
                 tmp2 * ( mat_( 0, 0 ) * mat_( 1, 1 ) - mat_( 1, 0 ) * mat_( 0, 1 ) );

      xComp *= real_c( 1 ) / jacDet_;
   }

   real_t evalDF( const Point3D& x, Matrix3r& DFx ) const override final
   {
      WALBERLA_UNUSED( x );
      DFx = mat_;
      return jacDet_;
   }

   bool verifyPointPairing( const Point3D&, const Point3D& ) const override final { return true; };

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override final
   {
      sendBuffer << Type::AFFINE_3D;
      for ( uint_t i = 0; i < 3; i++ )
      {
         for ( uint_t j = 0; j < 3; j++ )
         {
            sendBuffer << mat_( i, j );
         }
      }
      sendBuffer << vec_[0] << vec_[1] << vec_[2];
   }

   static void setMap( SetupPrimitiveStorage& setupStorage, const Matrix3r& mat, const Point3D& vec )
   {
      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(), std::make_shared< AffineMap3D >( mat, vec ) );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(), std::make_shared< AffineMap3D >( mat, vec ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(), std::make_shared< AffineMap3D >( mat, vec ) );
      }

      for ( auto it : setupStorage.getVertices() )
      {
         Vertex& vertex = *it.second;
         setupStorage.setGeometryMap( vertex.getID(), std::make_shared< AffineMap3D >( mat, vec ) );
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
      WALBERLA_ABORT( "AffineMap3D::evalDF unimplemented for 2D!" );
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const override final
   {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( DFinvx );
      WALBERLA_ABORT( "AffineMap3D::evalDFinv unimplemented for 2D!" );
   }
   ///@}

 private:
   /// matrix using in affine mapping
   Matrix3r mat_;

   /// translation vector
   Point3D vec_;

   /// value of Jacobian determinant
   real_t jacDet_;
};

} // end of namespace hyteg
