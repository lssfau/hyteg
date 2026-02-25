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
#pragma once

#include "core/DataTypes.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/communication/PackInfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFPackInfo.hpp"
#include "hyteg/p1functionspace/VertexDoFPackInfo.hpp"
#include "hyteg/types/Concepts.hpp"

namespace hyteg {

using walberla::uint_t;

template < typename ValueType >
class FunctionMemory;

template < typename DataType, typename PrimitiveType >
class PrimitiveDataID;

class Vertex;
class Edge;
class Face;
class Cell;
class PrimitiveStorage;
class PrimitiveID;

template < concepts::value_type ValueType >
class P2FunctionPackInfo : public communication::PackInfo
{
 public:
   P2FunctionPackInfo( uint_t                                                                  level,
                       std::array< PrimitiveDataID< FunctionMemory< ValueType >, Vertex >, 2 > dataIDsMacroVertex,
                       std::array< PrimitiveDataID< FunctionMemory< ValueType >, Edge >, 2 >   dataIDsMacroEdge,
                       std::array< PrimitiveDataID< FunctionMemory< ValueType >, Face >, 2 >   dataIDsMacroFace,
                       std::array< PrimitiveDataID< FunctionMemory< ValueType >, Cell >, 2 >   dataIDsMacroCell,
                       std::weak_ptr< PrimitiveStorage >                                       storage );

   void packVertexForEdge( const Vertex* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackEdgeFromVertex( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalVertexToEdge( const Vertex* sender, Edge* receiver ) const override;

   void packEdgeForVertex( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackVertexFromEdge( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalEdgeToVertex( const Edge* sender, Vertex* receiver ) const override;

   void packEdgeForFace( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackFaceFromEdge( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalEdgeToFace( const Edge* sender, Face* receiver ) const override;

   void packFaceForEdge( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackEdgeFromFace( Edge* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalFaceToEdge( const Face* sender, Edge* receiver ) const override;

   void packFaceForCell( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackCellFromFace( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalFaceToCell( const Face* sender, Cell* receiver ) const override;

   void packCellForFace( const Cell* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackFaceFromCell( Face* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalCellToFace( const Cell* sender, Face* receiver ) const override;

   void packVertexForCell( const Vertex* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackCellFromVertex( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalVertexToCell( const Vertex* sender, Cell* receiver ) const override;

   void packEdgeForCell( const Edge* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackCellFromEdge( Cell* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalEdgeToCell( const Edge* sender, Cell* receiver ) const override;

   void packFaceForVertex( const Face* sender, const PrimitiveID& receiver, walberla::mpi::SendBuffer& buffer ) const override;

   void unpackVertexFromFace( Vertex* receiver, const PrimitiveID& sender, walberla::mpi::RecvBuffer& buffer ) const override;

   void communicateLocalFaceToVertex( const Face* sender, Vertex* receiver ) const override;

 private:
   uint_t                                   level_;
   std::weak_ptr< hyteg::PrimitiveStorage > storage_;

   VertexDoFPackInfo< ValueType > vertexDoFPackInfo_;
   EdgeDoFPackInfo< ValueType >   edgeDoFPackInfo_;
};

} //namespace hyteg
