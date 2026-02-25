/*
 * Copyright (c) 2017-2026 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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
#include "P2FunctionPackInfo.hpp"

#include "hyteg/Algorithms.hpp"
#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/StencilDirections.hpp"
#include "hyteg/communication/DoFSpacePackInfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/indexing/DistanceCoordinateSystem.hpp"
#include "hyteg/indexing/LocalIDMappings.hpp"
#include "hyteg/memory/FunctionMemory.hpp"

namespace hyteg {

template < concepts::value_type ValueType >
P2FunctionPackInfo< ValueType >::P2FunctionPackInfo(
    uint_t                                                                  level,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Vertex >, 2 > dataIDsMacroVertex,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Edge >, 2 >   dataIDsMacroEdge,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Face >, 2 >   dataIDsMacroFace,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Cell >, 2 >   dataIDsMacroCell,
    std::weak_ptr< PrimitiveStorage >                                       storage )
: vertexDoFPackInfo_( level, dataIDsMacroVertex[0], dataIDsMacroEdge[0], dataIDsMacroFace[0], dataIDsMacroCell[0], storage )
, edgeDoFPackInfo_( level, dataIDsMacroVertex[1], dataIDsMacroEdge[1], dataIDsMacroFace[1], dataIDsMacroCell[1], storage )
{}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::packVertexForEdge( const Vertex*              sender,
                                                         const PrimitiveID&         receiver,
                                                         walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packVertexForEdge( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::unpackEdgeFromVertex( Edge*                      receiver,
                                                            const PrimitiveID&         sender,
                                                            walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackEdgeFromVertex( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalVertexToEdge( sender, receiver );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::packEdgeForVertex( const Edge*                sender,
                                                         const PrimitiveID&         receiver,
                                                         walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packEdgeForVertex( sender, receiver, buffer );
   edgeDoFPackInfo_.packEdgeForVertex( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::unpackVertexFromEdge( Vertex*                    receiver,
                                                            const PrimitiveID&         sender,
                                                            walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackVertexFromEdge( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackVertexFromEdge( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalEdgeToVertex( sender, receiver );
   edgeDoFPackInfo_.communicateLocalEdgeToVertex( sender, receiver );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::packEdgeForFace( const Edge*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packEdgeForFace( sender, receiver, buffer );
   edgeDoFPackInfo_.packEdgeForFace( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::unpackFaceFromEdge( Face*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackFaceFromEdge( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackFaceFromEdge( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalEdgeToFace( sender, receiver );
   edgeDoFPackInfo_.communicateLocalEdgeToFace( sender, receiver );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::packFaceForEdge( const Face*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packFaceForEdge( sender, receiver, buffer );
   edgeDoFPackInfo_.packFaceForEdge( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::unpackEdgeFromFace( Edge*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackEdgeFromFace( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackEdgeFromFace( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalFaceToEdge( sender, receiver );
   edgeDoFPackInfo_.communicateLocalFaceToEdge( sender, receiver );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::packFaceForCell( const Face*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packFaceForCell( sender, receiver, buffer );
   edgeDoFPackInfo_.packFaceForCell( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::unpackCellFromFace( Cell*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackCellFromFace( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackCellFromFace( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalFaceToCell( sender, receiver );
   edgeDoFPackInfo_.communicateLocalFaceToCell( sender, receiver );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::packCellForFace( const Cell*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packCellForFace( sender, receiver, buffer );
   edgeDoFPackInfo_.packCellForFace( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::unpackFaceFromCell( Face*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackFaceFromCell( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackFaceFromCell( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalCellToFace( sender, receiver );
   edgeDoFPackInfo_.communicateLocalCellToFace( sender, receiver );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::packVertexForCell( const Vertex*              sender,
                                                         const PrimitiveID&         receiver,
                                                         walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packVertexForCell( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::unpackCellFromVertex( Cell*                      receiver,
                                                            const PrimitiveID&         sender,
                                                            walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackCellFromVertex( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::communicateLocalVertexToCell( const Vertex* sender, Cell* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalVertexToCell( sender, receiver );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::packEdgeForCell( const Edge*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packEdgeForCell( sender, receiver, buffer );
   edgeDoFPackInfo_.packEdgeForCell( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::unpackCellFromEdge( Cell*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackCellFromEdge( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackCellFromEdge( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalEdgeToCell( sender, receiver );
   edgeDoFPackInfo_.communicateLocalEdgeToCell( sender, receiver );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::packFaceForVertex( const Face*                sender,
                                                         const PrimitiveID&         receiver,
                                                         walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packFaceForVertex( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::unpackVertexFromFace( Vertex*                    receiver,
                                                            const PrimitiveID&         sender,
                                                            walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackVertexFromFace( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P2FunctionPackInfo< ValueType >::communicateLocalFaceToVertex( const Face* sender, Vertex* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalFaceToVertex( sender, receiver );
}

template class P2FunctionPackInfo< double >;
template class P2FunctionPackInfo< float >;
template class P2FunctionPackInfo< int >;
template class P2FunctionPackInfo< long >;
template class P2FunctionPackInfo< long long >;

} // namespace hyteg
