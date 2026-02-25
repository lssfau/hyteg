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

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::P2FunctionPackInfoGeneric(
    uint_t                                                                  level,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Vertex >, 2 > dataIDsMacroVertex,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Edge >, 2 >   dataIDsMacroEdge,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Face >, 2 >   dataIDsMacroFace,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Cell >, 2 >   dataIDsMacroCell,
    std::weak_ptr< PrimitiveStorage >                                       storage )
: vertexDoFPackInfo_( level, dataIDsMacroVertex[0], dataIDsMacroEdge[0], dataIDsMacroFace[0], dataIDsMacroCell[0], storage )
, edgeDoFPackInfo_( level, dataIDsMacroVertex[1], dataIDsMacroEdge[1], dataIDsMacroFace[1], dataIDsMacroCell[1], storage )
{}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packVertexForEdge(
    const Vertex*              sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packVertexForEdge( sender, receiver, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackEdgeFromVertex(
    Edge*                      receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackEdgeFromVertex( receiver, sender, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalVertexToEdge(
    const Vertex* sender,
    Edge*         receiver ) const
{
   vertexDoFPackInfo_.communicateLocalVertexToEdge( sender, receiver );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packEdgeForVertex(
    const Edge*                sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packEdgeForVertex( sender, receiver, buffer );
   edgeDoFPackInfo_.packEdgeForVertex( sender, receiver, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackVertexFromEdge(
    Vertex*                    receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackVertexFromEdge( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackVertexFromEdge( receiver, sender, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalEdgeToVertex(
    const Edge* sender,
    Vertex*     receiver ) const
{
   vertexDoFPackInfo_.communicateLocalEdgeToVertex( sender, receiver );
   edgeDoFPackInfo_.communicateLocalEdgeToVertex( sender, receiver );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packEdgeForFace(
    const Edge*                sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packEdgeForFace( sender, receiver, buffer );
   edgeDoFPackInfo_.packEdgeForFace( sender, receiver, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackFaceFromEdge(
    Face*                      receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackFaceFromEdge( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackFaceFromEdge( receiver, sender, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalEdgeToFace(
    const Edge* sender,
    Face*       receiver ) const
{
   vertexDoFPackInfo_.communicateLocalEdgeToFace( sender, receiver );
   edgeDoFPackInfo_.communicateLocalEdgeToFace( sender, receiver );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packFaceForEdge(
    const Face*                sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packFaceForEdge( sender, receiver, buffer );
   edgeDoFPackInfo_.packFaceForEdge( sender, receiver, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackEdgeFromFace(
    Edge*                      receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackEdgeFromFace( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackEdgeFromFace( receiver, sender, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalFaceToEdge(
    const Face* sender,
    Edge*       receiver ) const
{
   vertexDoFPackInfo_.communicateLocalFaceToEdge( sender, receiver );
   edgeDoFPackInfo_.communicateLocalFaceToEdge( sender, receiver );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packFaceForCell(
    const Face*                sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packFaceForCell( sender, receiver, buffer );
   edgeDoFPackInfo_.packFaceForCell( sender, receiver, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackCellFromFace(
    Cell*                      receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackCellFromFace( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackCellFromFace( receiver, sender, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalFaceToCell(
    const Face* sender,
    Cell*       receiver ) const
{
   vertexDoFPackInfo_.communicateLocalFaceToCell( sender, receiver );
   edgeDoFPackInfo_.communicateLocalFaceToCell( sender, receiver );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packCellForFace(
    const Cell*                sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packCellForFace( sender, receiver, buffer );
   edgeDoFPackInfo_.packCellForFace( sender, receiver, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackFaceFromCell(
    Face*                      receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackFaceFromCell( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackFaceFromCell( receiver, sender, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalCellToFace(
    const Cell* sender,
    Face*       receiver ) const
{
   vertexDoFPackInfo_.communicateLocalCellToFace( sender, receiver );
   edgeDoFPackInfo_.communicateLocalCellToFace( sender, receiver );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packVertexForCell(
    const Vertex*              sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packVertexForCell( sender, receiver, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackCellFromVertex(
    Cell*                      receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackCellFromVertex( receiver, sender, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalVertexToCell(
    const Vertex* sender,
    Cell*         receiver ) const
{
   vertexDoFPackInfo_.communicateLocalVertexToCell( sender, receiver );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packEdgeForCell(
    const Edge*                sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packEdgeForCell( sender, receiver, buffer );
   edgeDoFPackInfo_.packEdgeForCell( sender, receiver, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackCellFromEdge(
    Cell*                      receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackCellFromEdge( receiver, sender, buffer );
   edgeDoFPackInfo_.unpackCellFromEdge( receiver, sender, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalEdgeToCell(
    const Edge* sender,
    Cell*       receiver ) const
{
   vertexDoFPackInfo_.communicateLocalEdgeToCell( sender, receiver );
   edgeDoFPackInfo_.communicateLocalEdgeToCell( sender, receiver );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packFaceForVertex(
    const Face*                sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packFaceForVertex( sender, receiver, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackVertexFromFace(
    Vertex*                    receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackVertexFromFace( receiver, sender, buffer );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalFaceToVertex(
    const Face* sender,
    Vertex*     receiver ) const
{
   vertexDoFPackInfo_.communicateLocalFaceToVertex( sender, receiver );
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packCellForEdge(
    const Cell*                sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   if constexpr ( std::is_same_v< VertexDoFPackInfoType, VertexDoFAdditivePackInfo< ValueType > > &&
                  std::is_same_v< EdgeDoFPackInfoType, EdgeDoFAdditivePackInfo< ValueType > > )
   {
      vertexDoFPackInfo_.packCellForEdge( sender, receiver, buffer );
      edgeDoFPackInfo_.packCellForEdge( sender, receiver, buffer );
   }
   else
   {
      WALBERLA_ABORT( "P2FunctionPackInfoGeneric::packCellForEdge(): cannot serve your request!" );
   }
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackEdgeFromCell(
    Edge*                      receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   if constexpr ( std::is_same_v< VertexDoFPackInfoType, VertexDoFAdditivePackInfo< ValueType > > &&
                  std::is_same_v< EdgeDoFPackInfoType, EdgeDoFAdditivePackInfo< ValueType > > )
   {
      vertexDoFPackInfo_.unpackEdgeFromCell( receiver, sender, buffer );
      edgeDoFPackInfo_.unpackEdgeFromCell( receiver, sender, buffer );
   }
   else
   {
      WALBERLA_ABORT( "P2FunctionPackInfoGeneric::unpackEdgeForCell(): cannot serve your request!" );
   }
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalCellToEdge(
    const Cell* sender,
    Edge*       receiver ) const
{
   if constexpr ( std::is_same_v< VertexDoFPackInfoType, VertexDoFAdditivePackInfo< ValueType > > &&
                  std::is_same_v< EdgeDoFPackInfoType, EdgeDoFAdditivePackInfo< ValueType > > )
   {
      vertexDoFPackInfo_.communicateLocalCellToEdge( sender, receiver );
      edgeDoFPackInfo_.communicateLocalCellToEdge( sender, receiver );
   }
   else
   {
      WALBERLA_ABORT( "P2FunctionPackInfoGeneric::communicateLocalCellToEdge(): cannot serve your request!" );
   }
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::packCellForVertex(
    const Cell*                sender,
    const PrimitiveID&         receiver,
    walberla::mpi::SendBuffer& buffer ) const
{
   if constexpr ( std::is_same_v< VertexDoFPackInfoType, VertexDoFAdditivePackInfo< ValueType > > &&
                  std::is_same_v< EdgeDoFPackInfoType, EdgeDoFAdditivePackInfo< ValueType > > )
   {
      vertexDoFPackInfo_.packCellForVertex( sender, receiver, buffer );
      edgeDoFPackInfo_.packCellForVertex( sender, receiver, buffer );
   }
   else
   {
      WALBERLA_ABORT( "P2FunctionPackInfoGeneric::packCellForVertex(): cannot serve your request!" );
   }
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::unpackVertexFromCell(
    Vertex*                    receiver,
    const PrimitiveID&         sender,
    walberla::mpi::RecvBuffer& buffer ) const
{
   if constexpr ( std::is_same_v< VertexDoFPackInfoType, VertexDoFAdditivePackInfo< ValueType > > &&
                  std::is_same_v< EdgeDoFPackInfoType, EdgeDoFAdditivePackInfo< ValueType > > )
   {
      vertexDoFPackInfo_.unpackVertexFromCell( receiver, sender, buffer );
      edgeDoFPackInfo_.unpackVertexFromCell( receiver, sender, buffer );
   }
   else
   {
      WALBERLA_ABORT( "P2FunctionPackInfoGeneric::unpackVertexFromCell(): cannot serve your request!" );
   }
}

template < concepts::value_type ValueType, typename VertexDoFPackInfoType, typename EdgeDoFPackInfoType >
void P2FunctionPackInfoGeneric< ValueType, VertexDoFPackInfoType, EdgeDoFPackInfoType >::communicateLocalCellToVertex(
    const Cell* sender,
    Vertex*     receiver ) const
{
   if constexpr ( std::is_same_v< VertexDoFPackInfoType, VertexDoFAdditivePackInfo< ValueType > > &&
                  std::is_same_v< EdgeDoFPackInfoType, EdgeDoFAdditivePackInfo< ValueType > > )
   {
      vertexDoFPackInfo_.communicateLocalCellToVertex( sender, receiver );
      edgeDoFPackInfo_.communicateLocalCellToVertex( sender, receiver );
   }
   else
   {
      WALBERLA_ABORT( "P2FunctionPackInfoGeneric::communicateLocalCellToVertex(): cannot serve your request!" );
   }
}

// Explicit instantiations for P2FunctionPackInfo
// clang-format off
template class P2FunctionPackInfoGeneric< double   , VertexDoFPackInfo< double    >, EdgeDoFPackInfo< double    > >;
template class P2FunctionPackInfoGeneric< float    , VertexDoFPackInfo< float     >, EdgeDoFPackInfo< float     > >;
template class P2FunctionPackInfoGeneric< int      , VertexDoFPackInfo< int       >, EdgeDoFPackInfo< int       > >;
template class P2FunctionPackInfoGeneric< long     , VertexDoFPackInfo< long      >, EdgeDoFPackInfo< long      > >;
template class P2FunctionPackInfoGeneric< long long, VertexDoFPackInfo< long long >, EdgeDoFPackInfo< long long > >;
// clang-format on

// Explicit instantiations for P2FunctionPackInfo
// clang-format off
template class P2FunctionPackInfoGeneric< double   , VertexDoFAdditivePackInfo< double    >, EdgeDoFAdditivePackInfo< double    > >;
template class P2FunctionPackInfoGeneric< float    , VertexDoFAdditivePackInfo< float     >, EdgeDoFAdditivePackInfo< float     > >;
template class P2FunctionPackInfoGeneric< int      , VertexDoFAdditivePackInfo< int       >, EdgeDoFAdditivePackInfo< int       > >;
template class P2FunctionPackInfoGeneric< long     , VertexDoFAdditivePackInfo< long      >, EdgeDoFAdditivePackInfo< long      > >;
template class P2FunctionPackInfoGeneric< long long, VertexDoFAdditivePackInfo< long long >, EdgeDoFAdditivePackInfo< long long > >;
// clang-format on

} // namespace hyteg
