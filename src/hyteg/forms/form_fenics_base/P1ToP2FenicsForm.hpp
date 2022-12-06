/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/fenics/ufc_traits.hpp"
#include "hyteg/forms/P2Form.hpp"

// P1 to P2

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "hyteg/forms/form_fenics_generated/p1_to_p2_divt.h"
#include "hyteg/forms/form_fenics_generated/p1_to_p2_tet_divt_tet.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace hyteg {

using walberla::real_c;

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P1ToP2FenicsForm : public Form
{
 public:
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrixr< 6, 3 > localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
   }

   void integrateVertexToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrixr< 6, 3 > localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 5, 0 );
      out[1] = localStiffnessMatrix( 5, 1 );
      out[2] = localStiffnessMatrix( 5, 2 );
   }

   // Method invoked by P1Elements3D (we return row of element matrix for vertex 0
   // at the vertex columns), so actually this should be better called
   // integrateVertexToVertex, but this would imply changing the P1 stuff, too.
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      Matrixr< 10, 4 > elMat;
      computeLocalStiffnessMatrix( coords, elMat );
      uint_t rowIdx = fenics::P2DoFMap[0][0];
      out[0]        = elMat( rowIdx, fenics::P2DoFMap[0][0] );
      out[1]        = elMat( rowIdx, fenics::P2DoFMap[1][1] );
      out[2]        = elMat( rowIdx, fenics::P2DoFMap[2][2] );
      out[3]        = elMat( rowIdx, fenics::P2DoFMap[3][3] );
   }

   real_t integrate( const std::array< Point3D, 4 >&     coords,
                     const P2Form::dofPosByVertexPair3D& cntrPos,
                     const P2Form::dofPosByVertexPair3D& leafPos ) const
   {
      Matrixr< 10, 4 > elMat;
      computeLocalStiffnessMatrix( coords, elMat );
      WALBERLA_ASSERT_LESS( leafPos[0], 4 );
      WALBERLA_ASSERT_LESS( leafPos[1], 4 );
      uint_t rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
      uint_t colIdx = fenics::P2DoFMap[leafPos[0]][leafPos[1]];

      return real_c( elMat( rowIdx, colIdx ) );
   }

   std::vector< real_t > integrate( const std::array< Point3D, 4 >&                    coords,
                                    const P2Form::dofPosByVertexPair3D&                cntrPos,
                                    const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const
   {
      WALBERLA_ABORT( "Missing implementation in P1ToP2FenicsForm" );
   }

   /// \brief Compute local element matrix and return it
   ///
   /// The method computes the local element matrix for the triangle
   /// given by the vertex coordinates. It returns the complete matrix.
   ///
   /// \param coords  The coordinates of the three vertices of the triangle
   /// \param elMat   On return is filled with the matrix entries
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrixr< 6, 3 > & elMat ) const
   {
      computeLocalStiffnessMatrix( coords, elMat );
   }

   /// \brief Compute local element matrix and return it
   ///
   /// The method computes the local element matrix for the tetrahedron
   /// given by the vertex coordinates. It returns the complete matrix.
   ///
   /// \param coords  The coordinates of the four vertices of the tetrahedron
   /// \param elMat   On return is filled with the matrix entries
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrixr< 10, 4 > & elMat ) const
   {
      computeLocalStiffnessMatrix( coords, elMat );
   }

   bool assemble2D() const override { return !std::is_same< UFCOperator2D, hyteg::fenics::NoAssemble >::value; }

   bool assemble3D() const override { return !std::is_same< UFCOperator3D, hyteg::fenics::NoAssemble >::value; }

   bool assembly2DDefined() const override { return !std::is_same< UFCOperator2D, hyteg::fenics::UndefinedAssembly >::value; }

   bool assembly3DDefined() const override { return !std::is_same< UFCOperator3D, hyteg::fenics::UndefinedAssembly >::value; }

   inline void setGeometryMap( const std::shared_ptr< GeometryMap > map ) const { WALBERLA_UNUSED( map ); }

 private:
   void computeLocalStiffnessMatrix( const std::array< Point3D, 3 >& coords, Matrixr< 6, 3 >& localStiffnessMatrix ) const
   {
      real_t fenicsCoords[6];
      fenicsCoords[0] = coords[0][0];
      fenicsCoords[1] = coords[0][1];
      fenicsCoords[2] = coords[1][0];
      fenicsCoords[3] = coords[1][1];
      fenicsCoords[4] = coords[2][0];
      fenicsCoords[5] = coords[2][1];
      UFCOperator2D gen;
      gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
   }

   void computeLocalStiffnessMatrix( const std::array< Point3D, 4 >& coords, Matrixr< 10, 4 >& localStiffnessMatrix ) const
   {
      real_t fenicsCoords[12];
      for ( int node = 0; node < 4; ++node )
      {
         for ( int dim = 0; dim < 3; ++dim )
         {
            fenicsCoords[node * 3 + dim] = coords[node][dim];
         }
      }
      UFCOperator3D gen;
      gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
   }
};

} // namespace hyteg
