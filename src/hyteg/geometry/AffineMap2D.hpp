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

/// The class implements a generic affine mapping in 2D
///
/// The affine mapping is characterised by a matrix \f$M\f$ and a vector \f$v\f$
/// and defined as
/// \f[
/// x \mapsto M x + v
/// \f]
class AffineMap2D : public GeometryMap
{
 public:
   AffineMap2D( const Matrix2r& mat, const Point2D& vec )
   : mat_( mat )
   , vec_( vec )
   {
      // precompute Jacobian determinant
      real_t tmp1 = +mat_( 0, 0 ) * mat_( 1, 1 );
      real_t tmp2 = -mat_( 1, 0 ) * mat_( 0, 1 );
      jacDet_     = tmp1 + tmp2;
   }

   AffineMap2D( walberla::mpi::RecvBuffer& recvBuffer )
   {
      for ( uint_t i = 0; i < 2; i++ )
      {
         for ( uint_t j = 0; j < 2; j++ )
         {
            recvBuffer >> mat_( i, j );
         }
      }
      recvBuffer >> vec_[0];
      recvBuffer >> vec_[1];
   }

   void evalF( const Point3D& xold, Point3D& xnew ) const override final
   {
      xnew[0] = mat_( 0, 0 ) * xold[0] + mat_( 0, 1 ) * xold[1] + vec_[0];
      xnew[1] = mat_( 1, 0 ) * xold[0] + mat_( 1, 1 ) * xold[1] + vec_[1];
      xnew[2] = real_c( 0 );
   }

   void evalFinv( const Point3D& xPhys, Point3D& xComp ) const
   {
      real_t tmp0 = real_c( 1 ) / jacDet_;
      real_t tmp1 = tmp0 * ( -vec_[0] + xPhys[0] );
      real_t tmp2 = tmp0 * ( -vec_[1] + xPhys[1] );

      xComp[0] = -mat_( 0, 1 ) * tmp2 + mat_( 1, 1 ) * tmp1;
      xComp[1] = mat_( 0, 0 ) * tmp2 - mat_( 1, 0 ) * tmp1;
   }

   void serializeSubClass( walberla::mpi::SendBuffer& sendBuffer ) const override final
   {
      sendBuffer << Type::AFFINE_2D;
      for ( uint_t i = 0; i < 2; i++ )
      {
         for ( uint_t j = 0; j < 2; j++ )
         {
            sendBuffer << mat_( i, j );
         }
      }
      sendBuffer << vec_[0] << vec_[1];
   }

   static void setMap( SetupPrimitiveStorage& setupStorage, const Matrix2r& mat, const Point2D& vec )
   {
      for ( auto it : setupStorage.getCells() )
      {
         Cell& cell = *it.second;
         setupStorage.setGeometryMap( cell.getID(), std::make_shared< AffineMap2D >( mat, vec ) );
      }

      for ( auto it : setupStorage.getFaces() )
      {
         Face& face = *it.second;
         setupStorage.setGeometryMap( face.getID(), std::make_shared< AffineMap2D >( mat, vec ) );
      }

      for ( auto it : setupStorage.getEdges() )
      {
         Edge& edge = *it.second;
         setupStorage.setGeometryMap( edge.getID(), std::make_shared< AffineMap2D >( mat, vec ) );
      }

      for ( auto it : setupStorage.getVertices() )
      {
         Vertex& vertex = *it.second;
         setupStorage.setGeometryMap( vertex.getID(), std::make_shared< AffineMap2D >( mat, vec ) );
      }
   }

   void evalDF( const Point3D& x, Matrix2r& DFx ) const override final
   {
      WALBERLA_UNUSED( x );
      DFx( 0, 0 ) = mat_( 0, 0 );
      DFx( 0, 1 ) = mat_( 0, 1 );
      DFx( 1, 0 ) = mat_( 1, 0 );
      DFx( 1, 1 ) = mat_( 1, 1 );
   }

   void evalDFinv( const Point3D& x, Matrix2r& DFinvx ) const override final
   {
      WALBERLA_UNUSED( x );
      DFinvx( 0, 0 ) = +mat_( 1, 1 );
      DFinvx( 0, 1 ) = -mat_( 0, 1 );
      DFinvx( 1, 0 ) = -mat_( 1, 0 );
      DFinvx( 1, 1 ) = +mat_( 0, 0 );
      DFinvx *= real_c( 1 ) / jacDet_;
   }

   bool verifyPointPairing( const Point3D& computationalCoordinates, const Point3D& physicalCoordinates ) const override final
   {
      return true;
   }

 private:
   /// matrix using in affine mapping
   Matrix2r mat_;

   /// translation vector
   Point2D vec_;

   /// value of Jacobian determinant
   real_t jacDet_;
};

} // end of namespace hyteg
