/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr.
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

#include "core/Abort.h"

#include "hyteg/forms/Form.hpp"

namespace hyteg {

class P2Form : public Form
{
 public:
   virtual ~P2Form() {}

   /// Describe a degree of freedom in a P2 element on a tetrahedron
   /// by two vertex indices. There are two cases:
   ///
   /// (a) Both vertex indices are identical, then this is the
   ///     index of the dof associated with this vertex.
   ///
   /// (b) The two vertex indices are different, then this is the
   ///     index of the dof associated with the midpoint of the
   ///     tet's edge given by those two vertices.
   typedef std::array< uint_t, 2 > dofPosByVertexPair3D;

   // ---------------------------
   //  2D versions for triangles
   // ---------------------------
   virtual void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   virtual void integrateEdgeToVertex( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   virtual void integrateVertexToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   virtual void integrateEdgeToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   virtual void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   virtual void integrateRow( const uint_t & row, const std::array< Point3D, 3 >& coords, Matrixr< 1, 6 >& elMat ) const;

   // ----------------------------
   //  3D versions for tetrahedra
   // ----------------------------
   virtual void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   virtual void integrateEdgeToVertex( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2LinearCombinationFenicsForm not implemented for 3D!" );
   }

   virtual void integrateVertexToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2LinearCombinationFenicsForm not implemented for 3D!" );
   }

   virtual void integrateEdgeToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2LinearCombinationFenicsForm not implemented for 3D!" );
   }

   virtual void integrateRow( const uint_t & row, const std::array< Point3D, 4 >& coords, Matrixr< 1, 10 >& elMat ) const;

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
   virtual real_t integrate( const std::array< Point3D, 4 >&     coords,
                             const P2Form::dofPosByVertexPair3D& cntrPos,
                             const P2Form::dofPosByVertexPair3D& leafPos ) const
   {
      return 0;  
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
   virtual std::vector< real_t > integrate( const std::array< Point3D, 4 >&                    coords,
                                            const P2Form::dofPosByVertexPair3D&                cntrPos,
                                            const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const
   {
      return std::vector< real_t >();  
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
   virtual void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const {}

 private:

   virtual void integrateRow0( const std::array< Point3D, 3 >& coords, Matrixr< 1, 6 >& elMat ) const;
   virtual void integrateRow0( const std::array< Point3D, 4 >& coords, Matrixr< 1, 10 >& elMat ) const;
};

} // namespace hyteg
