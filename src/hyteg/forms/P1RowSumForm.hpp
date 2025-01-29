/*
 * Copyright (c) 2020 Marcus Mohr.
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

#include <vector>

#include "core/debug/CheckFunctions.h"

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/forms/P1Form.hpp"

namespace hyteg {

/// \brief P1 form that writes the row sum of the underlying element matrix to the diagonal.
class P1RowSumForm : public P1Form
{
 public:
   P1RowSumForm()
   {
      WALBERLA_ABORT( "P1RowSumForm() only implemented for technical reasons. Please use P1RowSumForm( form ) instead!" );
   }

   explicit P1RowSumForm( const std::shared_ptr< P1Form >& form )
   : form_( form )
   {}

   ~P1RowSumForm() override = default;

   // ---------------------------
   //  2D versions for triangles
   // ---------------------------
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const override
   {
      out.setZero();
      Matrix3r elMat;
      form_->setGeometryMap( this->geometryMap_ );
      form_->integrateAll( coords, elMat );
      for ( uint_t j = 0; j < 3; j++ )
         out[0] += elMat( 0, j );
   }

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix3r& elMat ) const override
   {
      Matrix3r elementMatrix;
      form_->setGeometryMap( this->geometryMap_ );
      form_->integrateAll( coords, elementMatrix );
      elMat.setZero();
      for ( uint_t i = 0; i < 3; i++ )
      {
         real_t sum = 0;
         for ( uint_t j = 0; j < 3; j++ )
         {
            sum += elementMatrix( i, j );
         }
         elMat( i, i ) = sum;
      }
   }

   // ----------------------------
   //  3D versions for tetrahedra
   // ----------------------------
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const override
   {
      Matrix4r elementMatrix;
      integrateAll( coords, elementMatrix );
      out[0] = elementMatrix( 0, 0 );
      out[1] = elementMatrix( 0, 1 );
      out[2] = elementMatrix( 0, 2 );
      out[3] = elementMatrix( 0, 3 );
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix4r& elMat ) const override
   {
      Matrix4r elementMatrix;
      form_->setGeometryMap( geometryMap_ );
      form_->integrateAll( coords, elementMatrix );
      elMat.setZero();
      for ( uint_t i = 0; i < 4; i++ )
      {
         real_t sum = 0;
         for ( uint_t j = 0; j < 4; j++ )
         {
            sum += elementMatrix( i, j );
         }
         elMat( i, i ) = sum;
      }
   }

 private:

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrixr< 1, 3 >& elMat ) const override
   {
      Point3D row;
      integrate( coords, row );
      for ( int i = 0; i < 3; ++i )
      {
         elMat( 0, i ) = row[i];
      }
   }

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrixr< 1, 4 >& elMat ) const override
   {
      Point4D row;
      integrate( coords, row );
      for ( int i = 0; i < 4; ++i )
      {
         elMat( 0, i ) = row[i];
      }
   }

   std::shared_ptr< P1Form > form_;
};

} // namespace hyteg
