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

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

// P1
#include "hyteg/forms/form_fenics_generated/p1_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p1_div.h"
#include "hyteg/forms/form_fenics_generated/p1_div_K_grad.h"
#include "hyteg/forms/form_fenics_generated/p1_divt.h"
#include "hyteg/forms/form_fenics_generated/p1_mass.h"
#include "hyteg/forms/form_fenics_generated/p1_polar_laplacian.h"
#include "hyteg/forms/form_fenics_generated/p1_polar_mass.h"
#include "hyteg/forms/form_fenics_generated/p1_pspg.h"
#include "hyteg/forms/form_fenics_generated/p1_stokes_epsilon.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_div_tet.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_divt_tet.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_mass.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_pspg_tet.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_stokes_epsilon_tet.h"

// P2
#include "hyteg/forms/form_fenics_generated/p2_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_div.h"
#include "hyteg/forms/form_fenics_generated/p2_divt.h"
#include "hyteg/forms/form_fenics_generated/p2_mass.h"
#include "hyteg/forms/form_fenics_generated/p2_pspg.h"
#include "hyteg/forms/form_fenics_generated/p2_stokes_epsilon.h"
#include "hyteg/forms/form_fenics_generated/p2_stokes_full.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_div_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_divt_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_mass.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_pspg_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_stokes_epsilon_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_stokes_full_tet.h"

// P1 to P2
#include "hyteg/forms/form_fenics_generated/p1_to_p2_divt.h"
#include "hyteg/forms/form_fenics_generated/p1_to_p2_tet_divt_tet.h"

// P2 to P1
#include "hyteg/forms/form_fenics_generated/p2_to_p1_div.h"
#include "hyteg/forms/form_fenics_generated/p2_to_p1_tet_div_tet.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace hyteg {

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P2FenicsForm : public P2Form
{
 public:
   // ---------------------------
   //  2D versions for triangles
   // ---------------------------

   // Fenics P2 DoF ordering
   // 2         1---5--0
   // | \        \     |
   // |  \        \    |
   // 4   3        3   4
   // |    \        \  |
   // |     \        \ |
   // 0--5---1         2

   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const override
   {
      Matrix6r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
   }

   void integrateEdgeToVertex( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrix6r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 3 );
      out[1] = localStiffnessMatrix( 0, 4 );
      out[2] = localStiffnessMatrix( 0, 5 );
   }

   void integrateVertexToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrix6r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 5, 0 );
      out[1] = localStiffnessMatrix( 5, 1 );
      out[2] = localStiffnessMatrix( 5, 2 );
   }

   void integrateEdgeToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrix6r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 5, 3 );
      out[1] = localStiffnessMatrix( 5, 4 );
      out[2] = localStiffnessMatrix( 5, 5 );
   }

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const
   {
      computeLocalStiffnessMatrix( coords, elMat );
   }

   // ----------------------------
   //  3D versions for tetrahedra
   // ----------------------------
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const override
   {
      Matrix10r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
      out[3] = localStiffnessMatrix( 0, 3 );
   }

   void integrateEdgeToVertex( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2FenicsForm::integrateEdgeToVertex() not implemented for 3D!" );
   }

   void integrateVertexToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2FenicsForm::integrateVertexToEdge() not implemented for 3D!" );
   }

   void integrateEdgeToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2FenicsForm::integrateEdgeToEdge() not implemented for 3D!" );
   }

   // --------------------------------------------------
   //  3D versions for tetrahedra (using new interface)
   // --------------------------------------------------

   /// \brief Compute local element matrix and return a single entry
   ///
   /// The method computes the local element matrix for the tetrahedron
   /// given by the vertex coordinates. A single entry of this matrix
   /// is returned. Its row index corresponds to the dof given by cntrPos,
   /// the column index to that for the dof given by leafPos.
   ///
   /// \param coords   The coordinates of the four vertices of the tetrahedron
   /// \param cntrPos  degree of freedom specifying row index of entry
   /// \param leafPos  degree of freedom specifying column index of entry
   ///
   /// \return a single entry of the local element matrix
   real_t integrate( const std::array< Point3D, 4 >&     coords,
                     const P2Form::dofPosByVertexPair3D& cntrPos,
                     const P2Form::dofPosByVertexPair3D& leafPos ) const
   {
      Matrix10r elMat;
      computeLocalStiffnessMatrix( coords, elMat );
      uint_t rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
      uint_t colIdx = fenics::P2DoFMap[leafPos[0]][leafPos[1]];

      return walberla::real_c( elMat( rowIdx, colIdx ) );
   }

   /// \brief Compute local element matrix and return selected entries from a row
   ///
   /// The method computes the local element matrix for the tetrahedron
   /// given by the vertex coordinates. It returns selected entries from
   /// single row of this matrix. The row index corresponds to the dof
   /// given by cntrPos, the column indices in the row correspond to the
   /// dofs given by the leafPos array.
   ///
   /// \param coords   The coordinates of the four vertices of the tetrahedron
   /// \param cntrPos  degree of freedom specifying row index of entry
   /// \param leafPos  degrees of freedom specifying column indices of entries
   ///
   /// \return an array with entries of the local element matrix
   std::vector< real_t > integrate( const std::array< Point3D, 4 >&                         coords,
                                    const P2Form::dofPosByVertexPair3D&                     cntrPos,
                                    const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const
   {
      Matrix10r elMat;
      computeLocalStiffnessMatrix( coords, elMat );
      std::vector< real_t > matRow( leafPos.size(), 0 );

      uint_t rowIdx = fenics::P2DoFMap[cntrPos[0]][cntrPos[1]];
      uint_t colIdx = 0;

      for ( uint_t k = 0; k < leafPos.size(); ++k )
      {
         colIdx    = fenics::P2DoFMap[leafPos[k][0]][leafPos[k][1]];
         matRow[k] = walberla::real_c( elMat( rowIdx, colIdx ) );
      }

      return matRow;
   }

   /// \brief Compute local element matrix and return it
   ///
   /// The method computes the local element matrix for the tetrahedron
   /// given by the vertex coordinates. It returns the complete matrix.
   ///
   /// \Note Row and column indices correspond to FEniCS ordering for P2 element
   ///
   /// \param coords  The coordinates of the four vertices of the tetrahedron
   /// \param elMat   On return is filled with the matrix entries
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const
   {
      computeLocalStiffnessMatrix( coords, elMat );
   }

   // -------------

   bool assemble2D() const override { return !std::is_same< UFCOperator2D, hyteg::fenics::NoAssemble >::value; }

   bool assemble3D() const override { return !std::is_same< UFCOperator3D, hyteg::fenics::NoAssemble >::value; }

   bool assembly2DDefined() const override { return !std::is_same< UFCOperator2D, hyteg::fenics::UndefinedAssembly >::value; }

   bool assembly3DDefined() const override { return !std::is_same< UFCOperator3D, hyteg::fenics::UndefinedAssembly >::value; }

   inline void setGeometryMap( const std::shared_ptr< GeometryMap > map ) const { WALBERLA_UNUSED( map ); }

 private:
   void computeLocalStiffnessMatrix( const std::array< Point3D, 3 >& coords, Matrix6r& localStiffnessMatrix ) const
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

   void computeLocalStiffnessMatrix( const std::array< Point3D, 4 >& coords, Matrix10r& localStiffnessMatrix ) const
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
