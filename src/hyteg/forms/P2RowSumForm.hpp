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

#include "hyteg/forms/P2Form.hpp"

namespace hyteg {

/// \brief P2 form that writes the row sum of the underlying stiffness matrix to the diagonal.
class P2RowSumForm : public P2Form
{
 public:
   P2RowSumForm() = default;

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
      form_->integrateAll( coords, elMat );
      for ( uint_t j = 0; j < 6; j++ )
         out[2] += elMat( 5, j );
   }

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const
   {
      form_->integrateAll( coords, elMat );
      for ( uint_t i = 0; i < 6; i++ )
      {
         real_t sum = 0;
         for ( uint_t j = 0; j < 6; j++ )
         {
            sum += elMat( i, j );
         }
         elMat( i, i ) = sum;
      }
   }

   // ----------------------------
   //  3D versions for tetrahedra
   // ----------------------------
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2RowSum form not implemented for 3D." );
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
      WALBERLA_ABORT( "P2RowSum form not implemented for 3D." );
   }

   std::vector< real_t > integrate( const std::array< Point3D, 4 >&                    coords,
                                    const P2Form::dofPosByVertexPair3D&                cntrPos,
                                    const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const
   {
      WALBERLA_ABORT( "P2RowSum form not implemented for 3D." );
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const
   {
      WALBERLA_ABORT( "P2RowSum form not implemented for 3D." );
   }

   bool assemble2D() const override { return true; }

   bool assemble3D() const override { return false; }

   bool assembly2DDefined() const override { return true; }

   bool assembly3DDefined() const override { return false; }

 private:
   std::shared_ptr< P2Form > form_;
};

} // namespace hyteg
