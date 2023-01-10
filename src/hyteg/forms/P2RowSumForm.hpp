/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "hyteg/forms/P2Form.hpp"

namespace hyteg {

/// \brief P2 form that writes the row sum of the underlying stiffness matrix to the diagonal.
class P2RowSumForm : public P2Form
{
 public:
   P2RowSumForm()
   {
      WALBERLA_ABORT( "P2RowSumForm() only implemented for technical reasons. Please use P2RowSumForm( form ) instead!" );
   }

   P2RowSumForm( const std::shared_ptr< P2Form >& form )
   : form_( form )
   {}

   virtual ~P2RowSumForm() = default;

   // ---------------------------
   //  2D versions for triangles
   // ---------------------------
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      out.setAll( 0 );
      Matrix6r elMat;
      form_->setGeometryMap( geometryMap_ );
      form_->integrateAll( coords, elMat );
      for ( uint_t j = 0; j < 6; j++ )
         out[0] += elMat( 0, j );
   }

   void integrateEdgeToVertex( const std::array< Point3D, 3 >& coords, Point3D& out ) const { out.setAll( 0 ); }

   void integrateVertexToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const { out.setAll( 0 ); }

   void integrateEdgeToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      out.setAll( 0 );
      Matrix6r elMat;
      form_->setGeometryMap( geometryMap_ );
      form_->integrateAll( coords, elMat );
      for ( uint_t j = 0; j < 6; j++ )
         out[2] += elMat( 5, j );
   }

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const
   {
      Matrix6r stiffness;
      form_->setGeometryMap( geometryMap_ );
      form_->integrateAll( coords, stiffness );
      elMat.setAll( 0 );
      for ( uint_t i = 0; i < 6; i++ )
      {
         real_t sum = 0;
         for ( uint_t j = 0; j < 6; j++ )
         {
            sum += stiffness( i, j );
         }
         elMat( i, i ) = sum;
      }
   }

   // ----------------------------
   //  3D versions for tetrahedra
   // ----------------------------
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      Matrix10r localStiffnessMatrix;
      integrateAll( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
      out[3] = localStiffnessMatrix( 0, 3 );
   }

   void integrateEdgeToVertex( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2RowSum form not implemented for 3D." );
   }

   void integrateVertexToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2RowSum form not implemented for 3D." );
   }

   void integrateEdgeToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2RowSum form not implemented for 3D." );
   }

   // --------------------------------------------------
   //  3D versions for tetrahedra (using new interface)
   // --------------------------------------------------

   real_t integrate( const std::array< Point3D, 4 >&     coords,
                     const P2Form::dofPosByVertexPair3D& cntrPos,
                     const P2Form::dofPosByVertexPair3D& leafPos ) const
   {
      Matrix10r elMat;
      integrateAll( coords, elMat );
      int rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
      int colIdx = fenics::P2DoFMap[leafPos[0]][leafPos[1]];

      return walberla::real_c( elMat( rowIdx, colIdx ) );
   }

   std::vector< real_t > integrate( const std::array< Point3D, 4 >&                    coords,
                                    const P2Form::dofPosByVertexPair3D&                cntrPos,
                                    const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const
   {
      Matrix10r elMat;
      integrateAll( coords, elMat );
      std::vector< real_t > matRow( leafPos.size(), 0 );

      int rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
      int colIdx = 0;

      for ( uint_t k = 0; k < leafPos.size(); ++k )
      {
         colIdx    = fenics::P2DoFMap[leafPos[k][0]][leafPos[k][1]];
         matRow[k] = walberla::real_c( elMat( rowIdx, colIdx ) );
      }

      return matRow;
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const
   {
      Matrix10r stiffness;
      form_->setGeometryMap( geometryMap_ );
      form_->integrateAll( coords, stiffness );
      elMat.setAll( 0 );
      for ( uint_t i = 0; i < 10; i++ )
      {
         real_t sum = 0;
         for ( uint_t j = 0; j < 10; j++ )
         {
            sum += stiffness( i, j );
         }
         elMat( i, i ) = sum;
      }
   }

 private:
   std::shared_ptr< P2Form > form_;
};

} // namespace hyteg
