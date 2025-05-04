/*
 * Copyright (c) 2024 Andreas Burkhart.
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

#include <map>

#include "core/Environment.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitives/Vertex.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

enum MicroDoFType : size_t
{
   MacroInner = 0,

   MacroVertex0 = 1,
   MacroVertex1 = 2,
   MacroVertex2 = 3,
   MacroVertex3 = 4,

   MacroEdge0 = 5,
   MacroEdge1 = 6,
   MacroEdge2 = 7,
   MacroEdge3 = 8,
   MacroEdge4 = 9,
   MacroEdge5 = 10,

   MacroFace0 = 11,
   MacroFace1 = 12,
   MacroFace2 = 13,
   MacroFace3 = 14,
};

inline std::ostream& operator<<( std::ostream& os, const MicroDoFType type )
{
   switch ( type )
   {
   case MacroVertex0:
      return os << "MacroVertex0";
   case MacroVertex1:
      return os << "MacroVertex1";
   case MacroVertex2:
      return os << "MacroVertex2";
   case MacroVertex3:
      return os << "MacroVertex3";
   case MacroEdge0:
      return os << "MacroEdge0";
   case MacroEdge1:
      return os << "MacroEdge1";
   case MacroEdge2:
      return os << "MacroEdge2";
   case MacroEdge3:
      return os << "MacroEdge3";
   case MacroEdge4:
      return os << "MacroEdge4";
   case MacroEdge5:
      return os << "MacroEdge5";
   case MacroFace0:
      return os << "MacroFace0";
   case MacroFace1:
      return os << "MacroFace1";
   case MacroFace2:
      return os << "MacroFace2";
   case MacroFace3:
      return os << "MacroFace3";
   default:
      return os << "MacroInner";
   }
}

namespace celldof {

enum MicroCellBorderType : size_t
{
   BordersNoFace = 0,

   BordersFace0 = 1,
   BordersFace1 = 2,
   BordersFace2 = 4,
   BordersFace3 = 8,

   BordersFace0And1 = 1 + 2,
   BordersFace0And2 = 1 + 4,
   BordersFace0And3 = 1 + 8,
   BordersFace1And2 = 2 + 4,
   BordersFace1And3 = 2 + 8,
   BordersFace2And3 = 4 + 8,

   BordersFace0And1And2 = 1 + 2 + 4,
   BordersFace0And1And3 = 1 + 2 + 8,
   BordersFace0And2And3 = 1 + 4 + 8,
   BordersFace1And2And3 = 2 + 4 + 8,

   BordersAllFaces = 1 + 2 + 4 + 8
};

inline std::ostream& operator<<( std::ostream& os, const MicroCellBorderType type )
{
   switch ( type )
   {
   case BordersNoFace:
      return os << "BordersNoFace";
   case BordersFace0:
      return os << "BordersFace0";
   case BordersFace1:
      return os << "BordersFace1";
   case BordersFace2:
      return os << "BordersFace2";
   case BordersFace3:
      return os << "BordersFace3";
   case BordersFace0And1:
      return os << "BordersFace0And1";
   case BordersFace0And2:
      return os << "BordersFace0And2";
   case BordersFace0And3:
      return os << "BordersFace0And3";
   case BordersFace1And2:
      return os << "BordersFace1And2";
   case BordersFace1And3:
      return os << "BordersFace1And3";
   case BordersFace2And3:
      return os << "BordersFace2And3";
   case BordersFace0And1And2:
      return os << "BordersFace0And1And2";
   case BordersFace0And1And3:
      return os << "BordersFace0And1And3";
   case BordersFace0And2And3:
      return os << "BordersFace0And2And3";
   case BordersFace1And2And3:
      return os << "BordersFace1And2And3";
   default:
      return os << "BordersAllFaces";
   }
}

inline MicroCellBorderType operator|( MicroCellBorderType a, MicroCellBorderType b )
{
   return MicroCellBorderType( static_cast< size_t >( a ) | static_cast< size_t >( b ) );
}

inline MicroCellBorderType operator&( MicroCellBorderType a, MicroCellBorderType b )
{
   return MicroCellBorderType( static_cast< size_t >( a ) & static_cast< size_t >( b ) );
}

inline MicroCellBorderType operator^( MicroCellBorderType a, MicroCellBorderType b )
{
   return MicroCellBorderType( static_cast< size_t >( a ) ^ static_cast< size_t >( b ) );
}

MicroCellBorderType getMicroCellDoFTypeFromMicroCellBorderType( MicroDoFType type )
{
   switch ( type )
   {
   case MacroFace0:
      return BordersFace0;
      break;
   case MacroFace1:
      return BordersFace1;
      break;
   case MacroFace2:
      return BordersFace2;
      break;
   case MacroFace3:
      return BordersFace3;
      break;
   default:
      break;
   }

   return BordersNoFace;
}

// size_t is guaranteed to have 16 bits
enum MicroCellIndexType : size_t
{
   NoMatch = 0,

   Index0Equal0 = 1,
   Index1Equal0 = 2,
   Index2Equal0 = 4,

   Index0EqualFourthFaceIdentifier = 8,
   Index1EqualFourthFaceIdentifier = 16,
   Index2EqualFourthFaceIdentifier = 32,

   Index0And1Equal0 = 64,
   Index0And2Equal0 = 128,
   Index1And2Equal0 = 256,

   Index0And1EqualFourthFaceIdentifier = 512,
   Index0And2EqualFourthFaceIdentifier = 1024,
   Index1And2EqualFourthFaceIdentifier = 2048,

   IndexAllEqual0 = 4096,

   IndexAllEqualFourthFaceIdentifier = 8192,

   IndexBordersFace0And1 = 262,
   IndexBordersFace0And2 = 133,
   IndexBordersFace0And3 = 8708,
   IndexBordersFace1And2 = 67,
   IndexBordersFace1And3 = 9218,
   IndexBordersFace2And3 = 10241,

   IndexBordersFace0And1And2 = 4551,
   IndexBordersFace0And1And3 = 9998,
   IndexBordersFace0And2And3 = 10901,
   IndexBordersFace1And2And3 = 11363,

   IndexBordersAll = 16383
};

MicroCellBorderType getMicroCellBorderTypeFromIndexType( MicroCellIndexType indexType )
{
   switch ( indexType )
   {
   case Index2Equal0:
      return BordersFace0;
      break;
   case Index1Equal0:
      return BordersFace1;
      break;
   case Index0Equal0:
      return BordersFace2;
      break;
   case IndexAllEqualFourthFaceIdentifier:
      return BordersFace3;
      break;

   case IndexBordersFace0And1:
      return BordersFace0And1;
      break;
   case IndexBordersFace0And2:
      return BordersFace0And2;
      break;
   case IndexBordersFace0And3:
      return BordersFace0And3;
      break;
   case IndexBordersFace1And2:
      return BordersFace1And2;
      break;
   case IndexBordersFace1And3:
      return BordersFace1And3;
      break;
   case IndexBordersFace2And3:
      return BordersFace2And3;
      break;

   case IndexBordersFace0And1And2:
      return BordersFace0And1And2;
      break;
   case IndexBordersFace0And1And3:
      return BordersFace0And1And3;
      break;
   case IndexBordersFace0And2And3:
      return BordersFace0And2And3;
      break;
   case IndexBordersFace1And2And3:
      return BordersFace1And2And3;
      break;

   case IndexBordersAll:
      return BordersAllFaces;
      break;

   default:
      break;
   }

   return BordersNoFace;
}

std::map< MicroDoFType, hyteg::DoFType >
    getMicroDoFTypeToDoFTypeMapFromMacroCell( const std::shared_ptr< PrimitiveStorage >& storage,
                                              const Cell&                                cell,
                                              BoundaryCondition                          bc )
{
   // create MicroDoFType to DoFType map
   std::map< MicroDoFType, hyteg::DoFType > MicroDoFTypeToDoFTypeMap{ { MacroInner, hyteg::DoFType::Inner } };

   // get macro face DoFTypes
   std::vector< hyteg::PrimitiveID > neighborFaces;
   cell.getNeighborFaces( neighborFaces );

   for ( auto& faceId : neighborFaces )
   {
      Face& face = *( storage->getFace( faceId ) );

      auto localFaceID = cell.getLocalFaceID( faceId );

      hyteg::DoFType BCType = bc.getBoundaryType( face.getMeshBoundaryFlag() );

      MicroDoFTypeToDoFTypeMap[static_cast< MicroDoFType >( static_cast< size_t >( MacroFace0 ) +
                                                            static_cast< size_t >( localFaceID ) )] = BCType;
   }

   // get macro edge DoFTypes
   std::vector< hyteg::PrimitiveID > neighborEdges;
   cell.getNeighborEdges( neighborEdges );

   for ( auto& edgeId : neighborEdges )
   {
      Edge& edge = *( storage->getEdge( edgeId ) );

      auto localEdgeID = cell.getLocalEdgeID( edgeId );

      hyteg::DoFType BCType = bc.getBoundaryType( edge.getMeshBoundaryFlag() );

      MicroDoFTypeToDoFTypeMap[static_cast< MicroDoFType >( static_cast< size_t >( MacroEdge0 ) +
                                                            static_cast< size_t >( localEdgeID ) )] = BCType;
   }

   // get macro vertex DoFTypes
   std::vector< hyteg::PrimitiveID > neighborVertices;
   cell.getNeighborVertices( neighborVertices );

   for ( auto& vertexId : neighborVertices )
   {
      Vertex& vertex = *( storage->getVertex( vertexId ) );

      auto localVertexID = cell.getLocalVertexID( vertexId );

      hyteg::DoFType BCType = bc.getBoundaryType( vertex.getMeshBoundaryFlag() );

      MicroDoFTypeToDoFTypeMap[static_cast< MicroDoFType >( static_cast< size_t >( MacroVertex0 ) +
                                                            static_cast< size_t >( localVertexID ) )] = BCType;
   }

   return MicroDoFTypeToDoFTypeMap;
}

std::vector< MicroDoFType > getMicroDoFTypesFromMicroCell( const hyteg::celldof::CellType& cellType,
                                                           const hyteg::indexing::Index&   microCell,
                                                           uint_t                          level )
{
   // determine BorderType
   size_t indexType = static_cast< size_t >( NoMatch );

   // get fourth face index sum indicator 2^level - 1
   uint_t fourthFaceIndexSum = ( uint_t( 1 ) << level ) - 1;
   switch ( cellType )
   {
   case hyteg::celldof::CellType::BLUE_UP:
      fourthFaceIndexSum--;
      break;
   case hyteg::celldof::CellType::GREEN_UP:
      fourthFaceIndexSum--;
      break;
   case hyteg::celldof::CellType::WHITE_DOWN:
      fourthFaceIndexSum--;
      fourthFaceIndexSum--;
      break;
   case hyteg::celldof::CellType::BLUE_DOWN:
      fourthFaceIndexSum--;
      break;
   case hyteg::celldof::CellType::GREEN_DOWN:
      fourthFaceIndexSum--;
      break;
   default:
      break;
   }

   if ( microCell[0] == 0 )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index0Equal0 );
   }
   if ( microCell[1] == 0 )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index1Equal0 );
   }
   if ( microCell[2] == 0 )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index2Equal0 );
   }

   if ( microCell[0] == fourthFaceIndexSum )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index0EqualFourthFaceIdentifier );
   }
   if ( microCell[1] == fourthFaceIndexSum )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index1EqualFourthFaceIdentifier );
   }
   if ( microCell[2] == fourthFaceIndexSum )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index2EqualFourthFaceIdentifier );
   }

   if ( microCell[0] + microCell[1] == 0 )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index0And1Equal0 );
   }
   if ( microCell[0] + microCell[2] == 0 )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index0And2Equal0 );
   }
   if ( microCell[1] + microCell[2] == 0 )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index1And2Equal0 );
   }

   if ( microCell[0] + microCell[1] == fourthFaceIndexSum )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index0And1EqualFourthFaceIdentifier );
   }
   if ( microCell[0] + microCell[2] == fourthFaceIndexSum )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index0And2EqualFourthFaceIdentifier );
   }
   if ( microCell[1] + microCell[2] == fourthFaceIndexSum )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::Index1And2EqualFourthFaceIdentifier );
   }

   if ( microCell[0] + microCell[1] + microCell[2] == 0 )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::IndexAllEqual0 );
   }

   if ( microCell[0] + microCell[1] + microCell[2] == fourthFaceIndexSum )
   {
      indexType = indexType | static_cast< size_t >( MicroCellIndexType::IndexAllEqualFourthFaceIdentifier );
   }

   MicroCellBorderType borderType = getMicroCellBorderTypeFromIndexType( static_cast< MicroCellIndexType >( indexType ) );

   // determine MicroDofTypes
   // clang-format off
   switch ( cellType )
   {
   case hyteg::celldof::CellType::WHITE_UP:
      switch ( borderType )
      {
      case BordersFace0:
         return {MacroFace0, MacroFace0, MacroFace0, MacroInner, MacroFace0, MacroFace0, MacroFace0, MacroInner, MacroInner, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace1:
         return {MacroFace1, MacroFace1, MacroInner, MacroFace1, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace2:
         return {MacroFace2, MacroInner, MacroFace2, MacroFace2, MacroInner, MacroFace2, MacroInner, MacroFace2, MacroInner, MacroFace2, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner};
         break;
      case BordersFace3:
         return {MacroInner, MacroFace3, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroFace3, MacroInner, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner};
         break;   
      case BordersFace0And1:
         return {MacroEdge0, MacroEdge0, MacroFace0, MacroFace1, MacroEdge0, MacroFace0, MacroFace0, MacroFace1, MacroFace1, MacroInner, MacroFace0, MacroFace1, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersFace0And2:
         return {MacroEdge1, MacroFace0, MacroEdge1, MacroFace2, MacroFace0, MacroEdge1, MacroFace0, MacroFace2, MacroInner, MacroFace2, MacroFace0, MacroInner, MacroFace2, MacroInner, MacroInner};
         break; 
      case BordersFace0And3:
         return {MacroFace0, MacroEdge2, MacroEdge2, MacroFace3, MacroFace0, MacroFace0, MacroEdge2, MacroInner, MacroFace3, MacroFace3, MacroFace0, MacroInner, MacroInner, MacroFace3, MacroInner};
         break; 
      case BordersFace1And2:
         return {MacroEdge3, MacroFace1, MacroFace2, MacroEdge3, MacroFace1, MacroFace2, MacroInner, MacroEdge3, MacroFace1, MacroFace2, MacroInner, MacroFace1, MacroFace2, MacroInner, MacroInner};
         break; 
      case BordersFace1And3:
         return {MacroFace1, MacroEdge4, MacroFace3, MacroEdge4, MacroFace1, MacroInner, MacroFace3, MacroFace1, MacroEdge4, MacroFace3, MacroInner, MacroFace1, MacroInner, MacroFace3, MacroInner};
         break; 
      case BordersFace2And3:
         return {MacroFace2, MacroFace3, MacroEdge5, MacroEdge5, MacroInner, MacroFace2, MacroFace3, MacroFace2, MacroFace3, MacroEdge5, MacroInner, MacroInner, MacroFace2, MacroFace3, MacroInner};
         break;    
      case BordersFace0And1And2:
         return {MacroVertex0, MacroEdge0, MacroEdge1, MacroEdge3, MacroEdge0, MacroEdge1, MacroFace0, MacroEdge3, MacroFace1, MacroFace2, MacroFace0, MacroFace1, MacroFace2, MacroInner, MacroInner};
         break; 
      case BordersFace0And1And3:
         return {MacroEdge0, MacroVertex1, MacroEdge2, MacroEdge4, MacroEdge0, MacroFace0, MacroEdge2, MacroFace1, MacroEdge4, MacroFace3, MacroFace0, MacroFace1, MacroInner, MacroFace3, MacroInner};
         break; 
      case BordersFace0And2And3:
         return {MacroEdge1, MacroEdge2, MacroVertex2, MacroEdge5, MacroFace0, MacroEdge1, MacroEdge2, MacroFace2, MacroFace3, MacroEdge5, MacroFace0, MacroInner, MacroFace2, MacroFace3, MacroInner};
         break; 
      case BordersFace1And2And3:
         return {MacroEdge3, MacroEdge4, MacroEdge5, MacroVertex3, MacroFace1, MacroFace2, MacroFace3, MacroEdge3, MacroEdge4, MacroEdge5, MacroInner, MacroFace1, MacroFace2, MacroFace3, MacroInner};
         break;     
      case BordersAllFaces:
         return {MacroVertex0, MacroVertex1, MacroVertex2, MacroVertex3, MacroEdge0, MacroEdge1, MacroEdge2, MacroEdge3, MacroEdge4, MacroEdge5, MacroFace0, MacroFace1, MacroFace2, MacroFace3, MacroInner};
         break;                                                                          
      default:
         break;
      }
      break;
   case hyteg::celldof::CellType::BLUE_UP:
      switch ( borderType )
      {
      case BordersFace0:
         return {MacroFace0, MacroFace0, MacroFace0, MacroInner, MacroFace0, MacroFace0, MacroFace0, MacroInner, MacroInner, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace1:
         return {MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace2:
         return {MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace3:
         return {MacroInner, MacroFace3, MacroInner, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersFace0And1:
         return {MacroEdge0, MacroFace0, MacroFace0, MacroFace1, MacroFace0, MacroFace0, MacroFace0, MacroFace1, MacroInner, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersFace0And2:
         return {MacroFace0, MacroFace0, MacroEdge1, MacroInner, MacroFace0, MacroFace0, MacroFace0, MacroInner, MacroInner, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And3:
         return {MacroFace0, MacroEdge2, MacroFace0, MacroFace3, MacroFace0, MacroFace0, MacroFace0, MacroInner, MacroFace3, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And2:
         return {MacroFace1, MacroInner, MacroFace2, MacroFace1, MacroInner, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And3:
         return {MacroFace1, MacroFace3, MacroInner, MacroEdge4, MacroInner, MacroInner, MacroInner, MacroFace1, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace2And3:
         return {MacroInner, MacroFace3, MacroFace2, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;    
      case BordersFace0And1And2:
         return {MacroEdge0, MacroFace0, MacroEdge1, MacroFace1, MacroFace0, MacroFace0, MacroFace0, MacroFace1, MacroInner, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And1And3:
         return {MacroEdge0, MacroEdge2, MacroFace0, MacroEdge4, MacroFace0, MacroFace0, MacroFace0, MacroFace1, MacroFace3, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And2And3:
         return {MacroFace0, MacroEdge2, MacroEdge1, MacroFace3, MacroFace0, MacroFace0, MacroFace0, MacroInner, MacroFace3, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And2And3:
         return {MacroFace1, MacroFace3, MacroFace2, MacroEdge4, MacroInner, MacroInner, MacroInner, MacroFace1, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;    
      case BordersAllFaces:
         return {MacroEdge0, MacroEdge2, MacroEdge1, MacroEdge4, MacroFace0, MacroFace0, MacroFace0, MacroFace1, MacroFace3, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner};
         break;                                                                           
      default:
         break;
      }
      break;    
case hyteg::celldof::CellType::GREEN_UP:
      switch ( borderType )
      {
      case BordersFace0:
         return {MacroFace0, MacroFace0, MacroInner, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace1:
         return {MacroFace1, MacroInner, MacroFace1, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner};
         break;
      case BordersFace2:
         return {MacroInner, MacroFace2, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace3:
         return {MacroInner, MacroInner, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersFace0And1:
         return {MacroEdge0, MacroFace0, MacroFace1, MacroFace1, MacroFace0, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner};
         break;   
      case BordersFace0And2:
         return {MacroFace0, MacroEdge1, MacroInner, MacroFace2, MacroFace0, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And3:
         return {MacroFace0, MacroFace0, MacroFace3, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And2:
         return {MacroFace1, MacroFace2, MacroFace1, MacroEdge3, MacroInner, MacroFace1, MacroInner, MacroFace1, MacroFace2, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner};
         break; 
      case BordersFace1And3:
         return {MacroFace1, MacroInner, MacroEdge4, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner};
         break; 
      case BordersFace2And3:
         return {MacroInner, MacroFace2, MacroFace3, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;    
      case BordersFace0And1And2:
         return {MacroEdge0, MacroEdge1, MacroFace1, MacroEdge3, MacroFace0, MacroFace1, MacroInner, MacroFace1, MacroFace2, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner};
         break; 
      case BordersFace0And1And3:
         return {MacroEdge0, MacroFace0, MacroEdge4, MacroFace1, MacroFace0, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner};
         break; 
      case BordersFace0And2And3:
         return {MacroFace0, MacroEdge1, MacroFace3, MacroFace2, MacroFace0, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And2And3:
         return {MacroFace1, MacroFace2, MacroEdge4, MacroEdge3, MacroInner, MacroFace1, MacroInner, MacroFace1, MacroFace2, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner};
         break;     
      case BordersAllFaces:
         return {MacroEdge0, MacroEdge1, MacroEdge4, MacroEdge3, MacroFace0, MacroFace1, MacroInner, MacroFace1, MacroFace2, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner};
         break;                                                                             
      default:
         break;
      }
      break;      
case hyteg::celldof::CellType::WHITE_DOWN:
      switch ( borderType )
      {
      case BordersFace0:
         return {MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace1:
         return {MacroInner, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace2:
         return {MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace3:
         return {MacroInner, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersFace0And1:
         return {MacroFace0, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersFace0And2:
         return {MacroFace0, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And3:
         return {MacroFace0, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And2:
         return {MacroInner, MacroInner, MacroFace2, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And3:
         return {MacroInner, MacroFace3, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace2And3:
         return {MacroInner, MacroFace3, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;    
      case BordersFace0And1And2:
         return {MacroFace0, MacroInner, MacroFace2, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And1And3:
         return {MacroFace0, MacroFace3, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And2And3:
         return {MacroFace0, MacroFace3, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And2And3:
         return {MacroInner, MacroFace3, MacroFace2, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersAllFaces:
         // Not Possible
         break;                                                                               
      default:
         break;
      }
      break;      
case hyteg::celldof::CellType::BLUE_DOWN:
      switch ( borderType )
      {
      case BordersFace0:
         return {MacroInner, MacroInner, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace1:
         return {MacroFace1, MacroInner, MacroFace1, MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace2:
         return {MacroInner, MacroFace2, MacroFace2, MacroFace2, MacroInner, MacroInner, MacroFace2, MacroInner, MacroFace2, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner};
         break;
      case BordersFace3:
         return {MacroFace3, MacroFace3, MacroInner, MacroInner, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersFace0And1:
         return {MacroFace1, MacroInner, MacroFace1, MacroFace0, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersFace0And2:
         return {MacroInner, MacroFace2, MacroFace2, MacroEdge1, MacroInner, MacroInner, MacroFace2, MacroInner, MacroFace2, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner};
         break; 
      case BordersFace0And3:
         return {MacroFace3, MacroFace3, MacroInner, MacroFace0, MacroFace3, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And2:
         return {MacroFace1, MacroFace2, MacroEdge3, MacroFace2, MacroInner, MacroFace1, MacroFace2, MacroInner, MacroFace2, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner};
         break; 
      case BordersFace1And3:
         return {MacroEdge4, MacroFace3, MacroFace1, MacroInner, MacroFace3, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace2And3:
         return {MacroFace3, MacroEdge5, MacroFace2, MacroFace2, MacroFace3, MacroInner, MacroFace2, MacroInner, MacroFace2, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner};
         break;    
      case BordersFace0And1And2:
         return {MacroFace1, MacroFace2, MacroEdge3, MacroEdge1, MacroInner, MacroFace1, MacroFace2, MacroInner, MacroFace2, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner};
         break; 
      case BordersFace0And1And3:
         return {MacroEdge4, MacroFace3, MacroFace1, MacroFace0, MacroFace3, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And2And3:
         return {MacroFace3, MacroEdge5, MacroFace2, MacroEdge1, MacroFace3, MacroInner, MacroFace2, MacroInner, MacroFace2, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner};
         break; 
      case BordersFace1And2And3:
         return {MacroEdge4, MacroEdge5, MacroEdge3, MacroFace2, MacroFace3, MacroFace1, MacroFace2, MacroInner, MacroFace2, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner};
         break;               
      case BordersAllFaces:
         return {MacroEdge4, MacroEdge5, MacroEdge3, MacroEdge1, MacroFace3, MacroFace1, MacroFace2, MacroInner, MacroFace2, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner};
         break;                                                                   
      default:
         break;
      }
      break;      
case hyteg::celldof::CellType::GREEN_DOWN:
      switch ( borderType )
      {
      case BordersFace0:
         return {MacroFace0, MacroFace0, MacroInner, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace1:
         return {MacroInner, MacroInner, MacroFace1, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace2:
         return {MacroFace2, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;
      case BordersFace3:
         return {MacroInner, MacroFace3, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroFace3, MacroInner, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner};
         break;   
      case BordersFace0And1:
         return {MacroFace0, MacroFace0, MacroFace1, MacroInner, MacroFace0, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break;   
      case BordersFace0And2:
         return {MacroEdge1, MacroFace0, MacroInner, MacroFace2, MacroFace0, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And3:
         return {MacroFace0, MacroEdge2, MacroFace3, MacroFace3, MacroFace0, MacroInner, MacroFace3, MacroInner, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner};
         break; 
      case BordersFace1And2:
         return {MacroFace2, MacroInner, MacroFace1, MacroFace2, MacroInner, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace1And3:
         return {MacroInner, MacroFace3, MacroEdge4, MacroFace3, MacroInner, MacroInner, MacroFace3, MacroInner, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner};
         break; 
      case BordersFace2And3:
         return {MacroFace2, MacroFace3, MacroFace3, MacroEdge5, MacroInner, MacroInner, MacroFace3, MacroFace2, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner};
         break;    
      case BordersFace0And1And2:
         return {MacroEdge1, MacroFace0, MacroFace1, MacroFace2, MacroFace0, MacroInner, MacroInner, MacroFace2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
         break; 
      case BordersFace0And1And3:
         return {MacroFace0, MacroEdge2, MacroEdge4, MacroFace3, MacroFace0, MacroInner, MacroFace3, MacroInner, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner};
         break; 
      case BordersFace0And2And3:
         return {MacroEdge1, MacroEdge2, MacroFace3, MacroEdge5, MacroFace0, MacroInner, MacroFace3, MacroFace2, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner};
         break; 
      case BordersFace1And2And3:
         return {MacroFace2, MacroFace3, MacroEdge4, MacroEdge5, MacroInner, MacroInner, MacroFace3, MacroFace2, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner};
         break;      
      case BordersAllFaces:
         return {MacroEdge1, MacroEdge2, MacroEdge4, MacroEdge5, MacroFace0, MacroInner, MacroFace3, MacroFace2, MacroFace3, MacroFace3, MacroInner, MacroInner, MacroInner, MacroFace3, MacroInner};
         break;                                                                            
      default:
         break;
      }
      break;                                
   default:
      break;
   }
   
   // vertex0, vertex1, vertex2, vertex3, edge0, edge1, edge2, edge3, edge4, edge5, face0, face1, face2, face3, cell
   return {MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
   // clang-format on
}

std::vector< hyteg::DoFType >
    getP2DoFTypesFromMicroDoFTypesFEniCSOrdering3D( const std::map< MicroDoFType, hyteg::DoFType >& MicroDoFTypeToDoFTypeMap,
                                                    const std::vector< MicroDoFType >&              MicroDoFTypes )
{
   // MircoDofTypes ordering
   /// 3
   /// |\`\.
   /// | 9 `\.
   /// |  \   8
   /// 7  2 _  `\.
   /// |  /  `-6 `\.
   /// | 5      `\_`\
   /// 0------4------1
   // In list (dofs 0-9):
   //       0        1        2        3      4      5      6      7      8      9     10     11     12     13    14
   // vertex0, vertex1, vertex2, vertex3, edge0, edge1, edge2, edge3, edge4, edge5, face0, face1, face2, face3, cell
   /// FEniCS ordering:
   /// 3
   /// |\`\.
   /// | 4 `\.
   /// |  \   5
   /// 7  2 _  `\.
   /// |  /  `-6 `\.
   /// | 8      `\_`\
   /// 0------9------1
   // In list (dofs 0-9):
   //       0        1        2        3      4      5      6      7      8      9
   // vertex0, vertex1, vertex2, vertex3, edge5, edge4, edge2, edge3, edge1, edge0
   // cellDoFType is dropped as it is not needed for P2
   //
   // Hence we need to return 0,1,2,3,9,8,6,7,5,4 from MicroDoFTypes
   return { MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 0 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 1 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 2 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 3 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 9 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 8 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 6 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 7 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 5 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 4 ) ) };
}

std::vector< hyteg::DoFType > getP2DoFTypesFromMicroDoFTypesMicroDoFTypeOrdering3D(
    const std::map< MicroDoFType, hyteg::DoFType >& MicroDoFTypeToDoFTypeMap,
    const std::vector< MicroDoFType >&              MicroDoFTypes )
{
   // MircoDofTypes ordering
   /// 3
   /// |\`\.
   /// | 9 `\.
   /// |  \   8
   /// 7  2 _  `\.
   /// |  /  `-6 `\.
   /// | 5      `\_`\
   /// 0------4------1
   // In list (dofs 0-9):
   //       0        1        2        3      4      5      6      7      8      9     10     11     12     13    14
   // vertex0, vertex1, vertex2, vertex3, edge0, edge1, edge2, edge3, edge4, edge5, face0, face1, face2, face3, cell

   return { MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 0 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 1 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 2 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 3 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 4 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 5 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 6 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 7 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 8 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 9 ) ) };
}

} // namespace celldof

namespace facedof {

enum MicroFaceBorderType : size_t
{
   BordersNoEdge = 0,

   BordersEdge0 = 1,
   BordersEdge1 = 2,
   BordersEdge2 = 4,

   BordersEdge0And1 = 1 + 2,
   BordersEdge2And0 = 1 + 4,
   BordersEdge2And1 = 2 + 4,

   BordersAllEdges = 1 + 2 + 4
};

inline std::ostream& operator<<( std::ostream& os, const MicroFaceBorderType type )
{
   switch ( type )
   {
   case BordersNoEdge:
      return os << "BordersNoEdge";
   case BordersEdge0:
      return os << "BordersEdge0";
   case BordersEdge1:
      return os << "BordersEdge1";
   case BordersEdge2:
      return os << "BordersEdge2";
   case BordersEdge0And1:
      return os << "BordersEdge0And1";
   case BordersEdge2And0:
      return os << "BordersEdge2And0";
   case BordersEdge2And1:
      return os << "BordersEdge2And1";
   default:
      return os << "BordersAllEdges";
   }
}

inline MicroFaceBorderType operator|( MicroFaceBorderType a, MicroFaceBorderType b )
{
   return MicroFaceBorderType( static_cast< size_t >( a ) | static_cast< size_t >( b ) );
}

inline MicroFaceBorderType operator&( MicroFaceBorderType a, MicroFaceBorderType b )
{
   return MicroFaceBorderType( static_cast< size_t >( a ) & static_cast< size_t >( b ) );
}

inline MicroFaceBorderType operator^( MicroFaceBorderType a, MicroFaceBorderType b )
{
   return MicroFaceBorderType( static_cast< size_t >( a ) ^ static_cast< size_t >( b ) );
}

MicroFaceBorderType getMicroFaceBorderTypeFromMicroDoFType( MicroDoFType type )
{
   switch ( type )
   {
   case MacroEdge0:
      return BordersEdge0;
      break;
   case MacroEdge1:
      return BordersEdge1;
      break;
   case MacroEdge2:
      return BordersEdge2;
      break;

   default:
      break;
   }

   return BordersNoEdge;
}

std::map< MicroDoFType, hyteg::DoFType >
    getMicroDoFTypeToDoFTypeMapFromMacroFace( const std::shared_ptr< PrimitiveStorage >& storage,
                                              const Face&                                face,
                                              BoundaryCondition                          bc )
{
   // create MicroDoFType to DoFType map
   std::map< MicroDoFType, hyteg::DoFType > MicroDoFTypeToDoFTypeMap{ { MacroInner, hyteg::DoFType::Inner } };

   // get macro edge DoFTypes
   std::vector< hyteg::PrimitiveID > neighborEdges;
   face.getNeighborEdges( neighborEdges );

   for ( auto& edgeId : neighborEdges )
   {
      Edge& edge = *( storage->getEdge( edgeId ) );

      auto localEdgeID = face.edge_index( edgeId );

      hyteg::DoFType BCType = bc.getBoundaryType( edge.getMeshBoundaryFlag() );

      MicroDoFTypeToDoFTypeMap[static_cast< MicroDoFType >( static_cast< size_t >( MacroEdge0 ) +
                                                            static_cast< size_t >( localEdgeID ) )] = BCType;
   }

   // get macro vertex DoFTypes
   std::vector< hyteg::PrimitiveID > neighborVertices;
   face.getNeighborVertices( neighborVertices );

   for ( auto& vertexId : neighborVertices )
   {
      Vertex& vertex = *( storage->getVertex( vertexId ) );

      auto localVertexID = face.vertex_index( vertexId );

      hyteg::DoFType BCType = bc.getBoundaryType( vertex.getMeshBoundaryFlag() );

      MicroDoFTypeToDoFTypeMap[static_cast< MicroDoFType >( static_cast< size_t >( MacroVertex0 ) +
                                                            static_cast< size_t >( localVertexID ) )] = BCType;
   }

   return MicroDoFTypeToDoFTypeMap;
}

std::vector< MicroDoFType > getMicroDoFTypesFromMicroFace( const hyteg::facedof::FaceType& faceType,
                                                           const hyteg::indexing::Index&   microFace,
                                                           uint_t                          level )
{
   // get third edge index sum indicator 2^level - 1
   uint_t thirdEdgeIndexSum = ( uint_t( 1 ) << level ) - 1;

   switch ( faceType )
   {
   case facedof::FaceType::BLUE:
      thirdEdgeIndexSum--;
      break;
   default: // facedof::FaceType::GRAY
      break;
   }

   // get micro face border type
   MicroFaceBorderType borderType = BordersNoEdge;
   if ( microFace[1] == 0 )
   {
      borderType = borderType | BordersEdge0;
   }
   if ( microFace[0] == 0 )
   {
      borderType = borderType | BordersEdge1;
   }
   if ( microFace.sum() == thirdEdgeIndexSum )
   {
      borderType = borderType | BordersEdge2;
   }

   // clang-format off
   if ( borderType & BordersAllEdges )
   {
      switch ( faceType )
      {
      case facedof::FaceType::GRAY: {
         switch ( borderType )
         {
         case BordersEdge0:
            return {MacroEdge0, MacroEdge0, MacroInner, MacroEdge0, MacroInner, MacroInner, MacroInner};
            break;
         case BordersEdge1:
            return {MacroEdge1, MacroInner, MacroEdge1, MacroInner, MacroEdge1, MacroInner, MacroInner};
            break;
         case BordersEdge2:
            return {MacroInner, MacroEdge2, MacroEdge2, MacroInner, MacroInner, MacroEdge2, MacroInner};
            break;
         case BordersEdge0And1:
            return {MacroVertex0, MacroEdge0, MacroEdge1, MacroEdge0, MacroEdge1, MacroInner, MacroInner};
            break;
         case BordersEdge2And0:
            return {MacroEdge0, MacroVertex1, MacroEdge2, MacroEdge0, MacroInner, MacroEdge2, MacroInner};
            break;
         case BordersEdge2And1:
            return {MacroEdge1, MacroEdge2, MacroVertex2, MacroInner, MacroEdge1, MacroEdge2, MacroInner};
            break;
         case BordersAllEdges:
            return {MacroVertex0, MacroVertex1, MacroVertex2, MacroEdge0, MacroEdge1, MacroEdge2, MacroInner};
            break;
         default:
            break;
         }
      }
      break;

      case facedof::FaceType::BLUE: {
         switch ( borderType )
         {
         case BordersEdge0:
            return {MacroEdge0, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
            break;
         case BordersEdge1:
            return {MacroInner, MacroInner, MacroEdge1, MacroInner, MacroInner, MacroInner, MacroInner};
            break;
         case BordersEdge2:
            return {MacroInner, MacroEdge2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
            break;
         case BordersEdge0And1:
            return {MacroEdge0, MacroInner, MacroEdge1, MacroInner, MacroInner, MacroInner, MacroInner};
            break;
         case BordersEdge2And0:
            return {MacroEdge0, MacroEdge2, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
            break;
         case BordersEdge2And1:
            return {MacroInner, MacroEdge2, MacroEdge1, MacroInner, MacroInner, MacroInner, MacroInner};
            break;
         case BordersAllEdges:
            return {MacroEdge0, MacroEdge2, MacroEdge1, MacroInner, MacroInner, MacroInner, MacroInner};
            break;
         default:
            break;
         }
      }
      break;

      default:
         break;
      }
   }

   // vertex0, vertex1, vertex2, edge0, edge1, edge2, face
   return {MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner, MacroInner};
   // clang-format on
}

std::vector< hyteg::DoFType >
    getP2DoFTypesFromMicroDoFTypesFEniCSOrdering2D( const std::map< MicroDoFType, hyteg::DoFType >& MicroDoFTypeToDoFTypeMap,
                                                    const std::vector< MicroDoFType >&              MicroDoFTypes )
{
   // MircoDofTypes ordering
   ///    2
   ///    | \
   ///    |  \
   ///    |   \
   ///    4    5
   ///    |     \
   ///    |      \
   ///    0---3---1
   // In list (dofs 0-5):
   //       0        1        2      3      4      5     6
   // vertex0, vertex1, vertex2, edge0, edge1, edge2, face
   // P2 Fenics ordering
   ///    2
   ///    | \
   ///    |  \
   ///    |   \
   ///    4    3
   ///    |     \
   ///    |      \
   ///    0---5---1
   // In list (dofs 0-5):
   //       0        1        2      3      4      5
   // vertex0, vertex1, vertex2, edge2, edge1, edge0
   // faceDoFType is dropped as it is not needed for P2
   //
   // Hence we need to return 0,1,2,5,4,3 from MicroDoFTypes
   return { MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 0 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 1 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 2 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 5 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 4 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 3 ) ) };
}

std::vector< hyteg::DoFType > getP2DoFTypesFromMicroDoFTypesMicroDoFTypeOrdering2D(
    const std::map< MicroDoFType, hyteg::DoFType >& MicroDoFTypeToDoFTypeMap,
    const std::vector< MicroDoFType >&              MicroDoFTypes )
{
   // MircoDofTypes ordering
   ///    2
   ///    | \
   ///    |  \
   ///    |   \
   ///    4    5
   ///    |     \
   ///    |      \
   ///    0---3---1
   // In list (dofs 0-5):
   //       0        1        2      3      4      5     6
   // vertex0, vertex1, vertex2, edge0, edge1, edge2, face

   return { MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 0 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 1 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 2 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 3 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 4 ) ),
            MicroDoFTypeToDoFTypeMap.at( MicroDoFTypes.at( 5 ) ) };
}

} // namespace facedof
} // namespace hyteg
