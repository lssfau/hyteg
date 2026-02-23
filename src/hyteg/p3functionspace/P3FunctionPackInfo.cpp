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
#include "P3FunctionPackInfo.hpp"

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
P3FunctionPackInfo< ValueType >::P3FunctionPackInfo(
    uint_t                                                                  level,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Vertex >, 3 > dataIDsMacroVertex,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Edge >, 3 >   dataIDsMacroEdge,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Face >, 3 >   dataIDsMacroFace,
    std::array< PrimitiveDataID< FunctionMemory< ValueType >, Cell >, 3 >   dataIDsMacroCell,
    std::weak_ptr< PrimitiveStorage >                                       storage )
: vertexDoFPackInfo_( level, dataIDsMacroVertex[0], dataIDsMacroEdge[0], dataIDsMacroFace[0], dataIDsMacroCell[0], storage )
, edgeDoFPackInfoBlue_( level, dataIDsMacroVertex[1], dataIDsMacroEdge[1], dataIDsMacroFace[1], dataIDsMacroCell[1], storage )
, edgeDoFPackInfoRed_( level, dataIDsMacroVertex[2], dataIDsMacroEdge[2], dataIDsMacroFace[2], dataIDsMacroCell[2], storage )
{}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::packVertexForEdge( const Vertex*              sender,
                                                         const PrimitiveID&         receiver,
                                                         walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packVertexForEdge( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::unpackEdgeFromVertex( Edge*                      receiver,
                                                            const PrimitiveID&         sender,
                                                            walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackEdgeFromVertex( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalVertexToEdge( sender, receiver );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::packEdgeForVertex( const Edge*                sender,
                                                         const PrimitiveID&         receiver,
                                                         walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packEdgeForVertex( sender, receiver, buffer );

   uint_t vertexIdOnEdge = sender->vertex_index( receiver );

   // blue data (from receiver's perspective) is always the first we send
   if ( vertexIdOnEdge == 0 )
   {
      edgeDoFPackInfoBlue_.packEdgeForVertex( sender, receiver, buffer );
      edgeDoFPackInfoRed_.packEdgeForVertex( sender, receiver, buffer );
   }
   else
   {
      edgeDoFPackInfoRed_.packEdgeForVertex( sender, receiver, buffer );
      edgeDoFPackInfoBlue_.packEdgeForVertex( sender, receiver, buffer );
   }
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::unpackVertexFromEdge( Vertex*                    receiver,
                                                            const PrimitiveID&         sender,
                                                            walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackVertexFromEdge( receiver, sender, buffer );

   // blue data is always the first we receive
   edgeDoFPackInfoBlue_.unpackVertexFromEdge( receiver, sender, buffer );
   edgeDoFPackInfoRed_.unpackVertexFromEdge( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalEdgeToVertex( sender, receiver );

   uint_t vertexIdOnEdge = sender->vertex_index( receiver->getID() );

   // blue data (from receiver's perspective) is always the first we communicate
   if ( vertexIdOnEdge == 0 )
   {
      edgeDoFPackInfoBlue_.communicateLocalEdgeToVertex( sender, receiver );
      edgeDoFPackInfoRed_.communicateLocalEdgeToVertex( sender, receiver );
   }
   else
   {
      edgeDoFPackInfoBlue_.communicateLocalEdgeToVertex( sender, receiver );
      edgeDoFPackInfoRed_.communicateLocalEdgeToVertex( sender, receiver );
   }
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::packEdgeForFace( const Edge*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packEdgeForFace( sender, receiver, buffer );

   // we always send our blue data first and let face sort it out
   edgeDoFPackInfoBlue_.packEdgeForFace( sender, receiver, buffer );
   edgeDoFPackInfoRed_.packEdgeForFace( sender, receiver, buffer );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::unpackFaceFromEdge( Face*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackFaceFromEdge( receiver, sender, buffer );

   // frist data we received is sender's blue data, possibly need to adapt to receiver's layout
   uint_t edgeIndexOnFace = receiver->edge_index( sender );
   int    faceDir         = receiver->getEdgeOrientation()[edgeIndexOnFace];

   if ( faceDir == +1 )
   {
      edgeDoFPackInfoBlue_.unpackFaceFromEdge( receiver, sender, buffer );
      edgeDoFPackInfoRed_.unpackFaceFromEdge( receiver, sender, buffer );
   }
   else if ( faceDir == -1 )
   {
      edgeDoFPackInfoRed_.unpackFaceFromEdge( receiver, sender, buffer );
      edgeDoFPackInfoBlue_.unpackFaceFromEdge( receiver, sender, buffer );
   }
   else
   {
      WALBERLA_ABORT( "FaceBoundaryDirection = " << faceDir << "! Really?" );
   }
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalEdgeToFace( sender, receiver );

   // need to check marco-edge versus face's edge orientation
   uint_t edgeIndexOnFace = receiver->edge_index( sender->getID() );
   int    faceDir         = receiver->getEdgeOrientation()[edgeIndexOnFace];

   if ( faceDir == +1 )
   {
      edgeDoFPackInfoBlue_.communicateLocalEdgeToFace( sender, receiver );
      edgeDoFPackInfoRed_.communicateLocalEdgeToFace( sender, receiver );
   }
   else if ( faceDir == -1 )
   {
      edgeDoFPackInfoRed_.communicateLocalEdgeToFace( sender, receiver );
      edgeDoFPackInfoBlue_.communicateLocalEdgeToFace( sender, receiver );
   }
   else
   {
      WALBERLA_ABORT( "FaceBoundaryDirection = " << faceDir << "! Really?" );
   }
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::packFaceForEdge( const Face*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   vertexDoFPackInfo_.packFaceForEdge( sender, receiver, buffer );

   // we will send the edge's blue data first, so we need compare orientations
   // we send face's blue data first
   uint_t edgeIndexOnFace = sender->edge_index( receiver );
   int    faceDir         = sender->getEdgeOrientation()[edgeIndexOnFace];

   if ( faceDir == +1 ) // aligned
   {
      edgeDoFPackInfoBlue_.packFaceForEdge( sender, receiver, buffer );
      edgeDoFPackInfoRed_.packFaceForEdge( sender, receiver, buffer );
   }
   else if ( faceDir == -1 )
   {
      edgeDoFPackInfoRed_.packFaceForEdge( sender, receiver, buffer );
      edgeDoFPackInfoBlue_.packFaceForEdge( sender, receiver, buffer );
   }
   else
   {
      WALBERLA_ABORT( "FaceBoundaryDirection = " << faceDir << "! Really?" );
   }
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::unpackEdgeFromFace( Edge*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   vertexDoFPackInfo_.unpackEdgeFromFace( receiver, sender, buffer );

   // sender packed our blue data first
   edgeDoFPackInfoBlue_.unpackEdgeFromFace( receiver, sender, buffer );
   edgeDoFPackInfoRed_.unpackEdgeFromFace( receiver, sender, buffer );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const
{
   vertexDoFPackInfo_.communicateLocalFaceToEdge( sender, receiver );

   // need to check marco-edge versus face's edge orientation
   uint_t edgeIndexOnFace = sender->edge_index( receiver->getID() );
   int    faceDir         = sender->getEdgeOrientation()[edgeIndexOnFace];

   if ( faceDir == +1 ) // aligned
   {
      edgeDoFPackInfoBlue_.communicateLocalFaceToEdge( sender, receiver );
      edgeDoFPackInfoRed_.communicateLocalFaceToEdge( sender, receiver );
   }
   else if ( faceDir == -1 ) // swapped
   {
      edgeDoFPackInfoRed_.communicateLocalFaceToEdge( sender, receiver );
      edgeDoFPackInfoBlue_.communicateLocalFaceToEdge( sender, receiver );
   }
   else
   {
      WALBERLA_ABORT( "FaceBoundaryDirection = " << faceDir << "! Really?" );
   }
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::packFaceForCell( const Face*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template <>
void P3FunctionPackInfo< real_t >::packFaceForCell( const Face*                sender,
                                                    const PrimitiveID&         receiver,
                                                    walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::unpackCellFromFace( Cell*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template <>
void P3FunctionPackInfo< real_t >::unpackCellFromFace( Cell*                      receiver,
                                                       const PrimitiveID&         sender,
                                                       walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template <>
void P3FunctionPackInfo< real_t >::communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::packCellForFace( const Cell*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::unpackFaceFromCell( Face*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template <>
void P3FunctionPackInfo< real_t >::communicateLocalCellToFace( const Cell* sender, Face* receiver ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::packVertexForCell( const Vertex*, const PrimitiveID&, walberla::mpi::SendBuffer& ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::unpackCellFromVertex( Cell*, const PrimitiveID&, walberla::mpi::RecvBuffer& ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::communicateLocalVertexToCell( const Vertex*, Cell* ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::packEdgeForCell( const Edge*                sender,
                                                       const PrimitiveID&         receiver,
                                                       walberla::mpi::SendBuffer& buffer ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::unpackCellFromEdge( Cell*                      receiver,
                                                          const PrimitiveID&         sender,
                                                          walberla::mpi::RecvBuffer& buffer ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template < concepts::value_type ValueType >
void P3FunctionPackInfo< ValueType >::communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const
{
   WALBERLA_ABORT( "P3FunctionPackInfo does not support cells, yet!" );
}

template class P3FunctionPackInfo< double >;
template class P3FunctionPackInfo< float >;
template class P3FunctionPackInfo< int >;
template class P3FunctionPackInfo< long >;
template class P3FunctionPackInfo< long long >;

} // namespace hyteg
