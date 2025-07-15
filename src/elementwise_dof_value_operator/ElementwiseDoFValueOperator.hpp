/*
* Copyright (c) 2025 Andreas Burkhart
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

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {

void communicateDoFValue( const P1Function< real_t >&, const uint_t );
void communicateDoFValue( const P2Function< real_t >&, const uint_t );
void communicateDoFValue( const EdgeDoFFunction< real_t >&, const uint_t );


template < typename ValueType >
uint_t getDoFDataFromMicroFace( const P1Function< ValueType >& f,
                                const Face&                    face,
                                const indexing::Index&         microFaceIndex,
                                const facedof::FaceType&       faceType,
                                const uint_t                   level,
                                std::vector< ValueType >&      DoFValues,
                                uint_t                         offset )
{
   ValueType* VertexData = face.getData( f.getFaceDataID() )->getPointer( level );

   std::array< uint_t, 3 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroFace( microFaceIndex, faceType, level, vertexDoFIndices );

   // extract local data
   for ( uint_t k = 0; k < 3; k++ )
   {
      DoFValues[offset + k] = VertexData[vertexDoFIndices[k]];
   }

   return offset + 3;
}

template < typename ValueType >
uint_t getDoFDataFromMicroCell( const P1Function< ValueType >& f,
                                const Cell&                    cell,
                                const indexing::Index&         microCellIndex,
                                const celldof::CellType&       cellType,
                                const uint_t                   level,
                                std::vector< ValueType >&      DoFValues,
                                uint_t                         offset )
{
   ValueType* VertexData = cell.getData( f.getCellDataID() )->getPointer( level );

   std::array< uint_t, 4 > vertexDoFIndices;
   vertexdof::getVertexDoFDataIndicesFromMicroCell( microCellIndex, cellType, level, vertexDoFIndices );

   // extract local data
   for ( uint_t k = 0; k < 4; k++ )
   {
      DoFValues[offset + k] = VertexData[vertexDoFIndices[k]];
   }

   return offset + 4;
}

template < typename ValueType >
uint_t getDoFDataFromMicroFace( const EdgeDoFFunction< ValueType >& f,
                                const Face&                         face,
                                const indexing::Index&              microFaceIndex,
                                const facedof::FaceType&            faceType,
                                const uint_t                        level,
                                std::vector< ValueType >&           DoFValues,
                                uint_t                              offset )
{
   ValueType* EdgeData = face.getData( f.getFaceDataID() )->getPointer( level );

   std::array< uint_t, 3 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFaceIndex, faceType, level, edgeDoFIndices );

   // extract local data
   for ( uint_t k = 0; k < 3; k++ )
   {
      DoFValues[offset + k] = EdgeData[edgeDoFIndices[k]];
   }

   return offset + 3;
}

template < typename ValueType >
uint_t getDoFDataFromMicroCell( const EdgeDoFFunction< ValueType >& f,
                                const Cell&                         cell,
                                const indexing::Index&              microCellIndex,
                                const celldof::CellType&            cellType,
                                const uint_t                        level,
                                std::vector< ValueType >&           DoFValues,
                                uint_t                              offset )
{
   ValueType* EdgeData = cell.getData( f.getCellDataID() )->getPointer( level );

   std::array< uint_t, 6 > edgeDoFIndices;
   edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCellIndex, cellType, level, edgeDoFIndices );

   // extract local data
   for ( uint_t k = 0; k < 6; k++ )
   {
      DoFValues[offset + k] = EdgeData[edgeDoFIndices[k]];
   }

   return offset + 6;
}

template < typename ValueType >
uint_t getDoFDataFromMicroFace( const P2Function< ValueType >& f,
                                const Face&                    face,
                                const indexing::Index&         microFaceIndex,
                                const facedof::FaceType&       faceType,
                                const uint_t                   level,
                                std::vector< ValueType >&      DoFValues,
                                uint_t                         offset )
{
   offset = getDoFDataFromMicroFace( f.getVertexDoFFunction(), face, microFaceIndex, faceType, level, DoFValues, offset );
   return getDoFDataFromMicroFace( f.getEdgeDoFFunction(), face, microFaceIndex, faceType, level, DoFValues, offset );
}

template < typename ValueType >
uint_t getDoFDataFromMicroCell( const P2Function< ValueType >& f,
                                const Cell&                    cell,
                                const indexing::Index&         microCellIndex,
                                const celldof::CellType&       cellType,
                                const uint_t                   level,
                                std::vector< ValueType >&      DoFValues,
                                uint_t                         offset )
{
   offset = getDoFDataFromMicroCell( f.getVertexDoFFunction(), cell, microCellIndex, cellType, level, DoFValues, offset );
   return getDoFDataFromMicroCell( f.getEdgeDoFFunction(), cell, microCellIndex, cellType, level, DoFValues, offset );
}

// ########################################
// #### number of local DoFs on a face ####
// ########################################

// default case
template < class Dummy >
struct functionLocalDoFsFace
{
   static const size_t value = 0;
};

// P1Function
template <>
struct functionLocalDoFsFace< vertexdof::VertexDoFFunction< real_t > >
{
   static const size_t value = 3;
};

// EdgeDoFFunction
template <>
struct functionLocalDoFsFace< EdgeDoFFunction< real_t > >
{
   static const size_t value = 3;
};

// P2Function
template <>
struct functionLocalDoFsFace< P2Function< real_t > >
{
   static const size_t value = 6;
};

// ########################################
// #### number of local DoFs on a cell ####
// ########################################

// default case
template < class Dummy >
struct functionLocalDoFsCell
{
   static const size_t value = 0;
};

// P1Function
template <>
struct functionLocalDoFsCell< vertexdof::VertexDoFFunction< real_t > >
{
   static const size_t value = 4;
};

// EdgeDoFFunction
template <>
struct functionLocalDoFsCell< EdgeDoFFunction< real_t > >
{
   static const size_t value = 6;
};

// P2Function
template <>
struct functionLocalDoFsCell< P2Function< real_t > >
{
   static const size_t value = 10;
};

} // namespace hyteg