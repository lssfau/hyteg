/*
 * Copyright (c) 2025-2026 Marcus Mohr.
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

#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"
#include "hyteg/volumedofspace/VolumeDoFIndexing.hpp"

namespace hyteg {

// ==============
//  DG1Operators
// ==============

/// Generate an identity (sub-)matrix for a DG1Function
///
/// Note: the flag argument will be ignored
inline void saveIdentityOperator( const DG1Function< idx_t >&                 numerator,
                                  const std::shared_ptr< SparseMatrixProxy >& mat,
                                  size_t                                      level,
                                  DoFType                                     flag )
{
   WALBERLA_UNUSED( flag );

   const auto storage = numerator.getStorage();

   if ( storage->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Missing 3D Implementation in DG1Petsc.hpp!" );
   }

   else
   {
      for ( const auto& pid : storage->getFaceIDs() )
      {
         const auto dstPolyDegree = uint_c( numerator.getDGFunction()->polynomialDegree( pid ) );
         const auto numDstDofs    = numerator.getDGFunction()->basis()->numDoFsPerElement( 2, dstPolyDegree );
         auto       dstDofMemory  = numerator.getDGFunction()->volumeDoFFunction()->dofMemory( pid, level );
         const auto dstMemLayout  = numerator.getDGFunction()->volumeDoFFunction()->memoryLayout();

         WALBERLA_ASSERT_EQUAL( dstPolyDegree, 1 );
         WALBERLA_ASSERT_EQUAL( numDstDofs, 3 );

         for ( const auto faceType : facedof::allFaceTypes )
         {
            for ( auto itFace = facedof::macroface::Iterator( level, faceType ).begin();
                  itFace != facedof::macroface::Iterator( level, faceType ).end();
                  ++itFace )
            {
               indexing::Index elementIdx = *itFace;

               // This object does the heavy lifting of computing all required coordinates and normals.
               volumedofspace::indexing::ElementNeighborInfo neighborInfo = volumedofspace::indexing::ElementNeighborInfo(
                   elementIdx, faceType, level, numerator.getBoundaryCondition(), pid, storage );

               for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
               {
                  const auto globalDiagonalIdx = dstDofMemory[volumedofspace::indexing::index(
                      elementIdx.x(), elementIdx.y(), faceType, dstDofIdx, numDstDofs, level, dstMemLayout )];

                  mat->addValue( globalDiagonalIdx, globalDiagonalIdx, real_c( 1 ) );
               }
            }
         }
      }
   }
}

} // namespace hyteg
