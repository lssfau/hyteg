/*
 * Copyright (c) 2017-2025 Marcus Mohr, Nils Kohl
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

#include "hyteg/p1functionspace/P1Elements.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"
#include "hyteg/p1functionspace/globalIndices.hpp"

namespace hyteg {
namespace p1 {

using walberla::real_t;

/// compute product of element local vector with element matrix
///
/// \param level          level on which we operate in mesh hierarchy
/// \param microFace      index associated with the current element = micro-face
/// \param fType          type of micro-face (GRAY or BLUE)
/// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
/// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
/// \param elMat          the 3x3 element matrix to be multiplied
/// \param alpha          scaling factor that is applied to the local result vector
///
/// \note The src and dst data arrays must not be identical.
static inline void localMatrixVectorMultiply2D( const uint_t           level,
                                  const indexing::Index& microFace,
                                  facedof::FaceType      fType,
                                  const real_t* const    srcVertexData,
                                  real_t* const          dstVertexData,
                                  const Matrix3r&        elMat,
                                  const real_t&          alpha )
{
   WALBERLA_ASSERT_UNEQUAL( srcVertexData, dstVertexData );

   // obtain data indices of dofs associated with micro-face
   std::array< uint_t, 3 > vertexDoFIndices;
   p1::getGlobalIndices2D(fType, level, microFace, vertexDoFIndices);

   // assemble local element vector
   Point3D elVecOld, elVecNew;
   for ( int k = 0; k < 3; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[uint_c( k )]];
   }

   // apply matrix (operator locally)
   elVecNew = alpha * ( elMat * elVecOld );

   // redistribute result from "local" to "global vector"
   for ( int k = 0; k < 3; ++k )
   {
      dstVertexData[vertexDoFIndices[uint_c( k )]] += elVecNew[k];
   }
}

/// compute product of element local vector with element matrix
///
/// \param level          level on which we operate in mesh hierarchy
/// \param microCell      index associated with the current element = micro-cell
/// \param cType          type of micro-cell (WHITE_UP, BLUE_DOWN, ...)
/// \param srcVertexData  pointer to DoF data on micro-vertices (for reading data)
/// \param dstVertexData  pointer to DoF data on micro-vertices (for writing data)
/// \param elMat          the 4x4 element matrix to be multiplied
/// \param alpha          scaling factor that is applied to the local result vector
///
/// \note The src and dst data arrays must not be identical.
static inline void localMatrixVectorMultiply3D( const uint_t            level,
                                  const indexing::Index&  microCell,
                                  const celldof::CellType cType,
                                  const real_t* const     srcVertexData,
                                  real_t* const           dstVertexData,
                                  const Matrix4r&         elMat,
                                  const real_t&           alpha )
{
   // obtain data indices of dofs associated with micro-cell
   std::array< uint_t, 4 > vertexDoFIndices;
   // p1::getGlobalIndices3D(cType, level, microCell, vertexDoFIndices);

   // assemble local element vector
   Point4D elVecOld, elVecNew;
   for ( int k = 0; k < 4; ++k )
   {
      elVecOld[k] = srcVertexData[vertexDoFIndices[uint_c( k )]];
   }

   // local matvec
   // elVecNew = elMat * elVecOld;
   elVecNew[0] = elMat(0,0) * elVecOld[0] + elMat(0,1) * elVecOld[1] + elMat(0,2) * elVecOld[2] + elMat(0,3) * elVecOld[3];
   elVecNew[1] = elMat(1,0) * elVecOld[0] + elMat(1,1) * elVecOld[1] + elMat(1,2) * elVecOld[2] + elMat(1,3) * elVecOld[3];
   elVecNew[2] = elMat(2,0) * elVecOld[0] + elMat(2,1) * elVecOld[1] + elMat(2,2) * elVecOld[2] + elMat(2,3) * elVecOld[3];
   elVecNew[3] = elMat(3,0) * elVecOld[0] + elMat(3,1) * elVecOld[1] + elMat(3,2) * elVecOld[2] + elMat(3,3) * elVecOld[3];
   // elVecNew[0] = elMat(0,0) * elVecOld[0];
   // elVecNew[1] = elMat(1,1) * elVecOld[1];
   // elVecNew[2] = elMat(2,2) * elVecOld[2];
   // elVecNew[3] = elMat(3,3) * elVecOld[3];

   // redistribute result from "local" to "global vector"
   for ( int k = 0; k < 4; ++k )
   {
      dstVertexData[vertexDoFIndices[uint_c( k )]] += alpha * elVecNew[k];
   }
}

template < class P1Form >
static inline void assembleLocalElementMatrix2D( const Face&            face,
                                   uint_t                 level,
                                   const indexing::Index& microFace,
                                   facedof::FaceType      fType,
                                   P1Form                 form,
                                   Matrix3r&              elMat )
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 3 > verts = facedof::macroface::getMicroVerticesFromMicroFace( microFace, fType );
   std::array< Point3D, 3 >         coords;
   for ( uint_t k = 0; k < 3; ++k )
   {
      coords[k] = vertexdof::macroface::coordinateFromIndex( level, face, verts[k] );
   }

   // assemble local element matrix
   form.setGeometryMap( face.getGeometryMap() );
   form.integrateAll( coords, elMat );
}

template < class P1Form >
static inline void assembleLocalElementMatrix3D( const Cell&            cell,
                                   uint_t                 level,
                                   const indexing::Index& microCell,
                                   celldof::CellType      cType,
                                   P1Form                 form,
                                   Matrix4r&              elMat )
{
   // determine coordinates of vertices of micro-element
   std::array< indexing::Index, 4 > verts = celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cType );
   std::array< Point3D, 4 >         coords;
   for ( uint_t k = 0; k < 4; ++k )
   {
      coords[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
   }

   // assemble local element matrix
   form.setGeometryMap( cell.getGeometryMap() );
   form.integrateAll( coords, elMat );
}

} // namespace p1
} // namespace hyteg
