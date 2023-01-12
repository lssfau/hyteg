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

   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const override;

   void integrateEdgeToVertex( const std::array< Point3D, 3 >& coords, Point3D& out ) const override;

   void integrateVertexToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const override;

   void integrateEdgeToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const override;

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const override;

   // ----------------------------
   //  3D versions for tetrahedra
   // ----------------------------
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const override;

   void integrateEdgeToVertex( const std::array< Point3D, 4 >&, Point4D& ) const override
   {
      WALBERLA_ABORT( "P2FenicsForm::integrateEdgeToVertex() not implemented for 3D!" );
   }

   void integrateVertexToEdge( const std::array< Point3D, 4 >&, Point4D& ) const override
   {
      WALBERLA_ABORT( "P2FenicsForm::integrateVertexToEdge() not implemented for 3D!" );
   }

   void integrateEdgeToEdge( const std::array< Point3D, 4 >&, Point4D& ) const override
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
                     const P2Form::dofPosByVertexPair3D& leafPos ) const override;

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
   std::vector< real_t > integrate( const std::array< Point3D, 4 >&                    coords,
                                    const P2Form::dofPosByVertexPair3D&                cntrPos,
                                    const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const override;

   /// \brief Compute local element matrix and return it
   ///
   /// The method computes the local element matrix for the tetrahedron
   /// given by the vertex coordinates. It returns the complete matrix.
   ///
   /// \Note Row and column indices correspond to FEniCS ordering for P2 element
   ///
   /// \param coords  The coordinates of the four vertices of the tetrahedron
   /// \param elMat   On return is filled with the matrix entries
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const override;

   inline void setGeometryMap( const std::shared_ptr< GeometryMap > map ) const { WALBERLA_UNUSED( map ); }

 private:
   void computeLocalStiffnessMatrix( const std::array< Point3D, 3 >& coords, Matrix6r& localStiffnessMatrix ) const;

   void computeLocalStiffnessMatrix( const std::array< Point3D, 4 >& coords, Matrix10r& localStiffnessMatrix ) const;
};

} // namespace hyteg
