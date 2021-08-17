/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "hyteg/facedofspace/FaceDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"

namespace hyteg {
namespace facedof {
namespace macroedge {

template < typename ValueType >
inline void
    enumerate( const uint_t& Level, Edge& edge, const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId, uint_t& num )
{
   ValueType*                      dst = edge.getData( dstId )->getPointer( Level );
   std::vector< stencilDirection > dirs;
   dirs.push_back( stencilDirection::CELL_GRAY_SE );
   if ( edge.getNumHigherDimNeighbors() == 2 )
   {
      dirs.push_back( stencilDirection::CELL_GRAY_NE );
   }
   for ( auto dir : dirs )
   {
      for ( uint_t i = 1; i < ( hyteg::levelinfo::num_microvertices_per_edge( Level ) - 2 ); ++i )
      {
         dst[facedof::macroedge::indexFaceFromVertex( Level, i, dir )] = num;
         ++num;
      }
   }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                                              Level,
                         Edge&                                                                      edge,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                edgeMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >& srcMemoryIds,
                         const std::function< ValueType( const hyteg::Point3D&, const std::vector< ValueType >& f ) >& expr,
                         const std::shared_ptr< PrimitiveStorage >&                                                    storage )
{
   auto edgeMemory = edge.getData( edgeMemoryId )->getPointer( Level );

   uint_t rowsize = levelinfo::num_microvertices_per_edge( Level );

   Point3D x;
   Point3D xBlend;
   Point3D x0 = edge.getCoordinates()[0];
   Point3D dx = edge.getDirection() / ( walberla::real_c( rowsize - 1 ) );

   Face*  face0              = storage->getFace( edge.neighborFaces()[0].getID() );
   uint_t outerVertexOnFace0 = face0->vertex_index( face0->get_vertex_opposite_to_edge( edge.getID() ) );

   Point3D face0d0 = dx;
   Point3D face0d2 = ( face0->getCoordinates()[outerVertexOnFace0] - x0 ) / ( walberla::real_c( rowsize - 1 ) );

   std::vector< ValueType* > srcPtr;
   for ( auto src : srcMemoryIds )
   {
      srcPtr.push_back( edge.getData( src )->getPointer( Level ) );
   }

   std::vector< ValueType > srcVector( srcMemoryIds.size() );

   // gray south cells
   x = x0 + 1.0 / 3.0 * ( face0d0 + face0d2 ) + dx;
   for ( size_t i = 1; i < rowsize - 2; ++i )
   {
      for ( size_t k = 0; k < srcPtr.size(); ++k )
      {
         srcVector[k] = srcPtr[k][facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )];
      }

      face0->getGeometryMap()->evalF( x, xBlend );
      edgeMemory[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )] = expr( xBlend, srcVector );
      x += dx;
   }

   if ( edge.getNumNeighborFaces() == 2 )
   {
      // gray north cells
      Face*  face1              = storage->getFace( edge.neighborFaces()[1].getID() );
      uint_t outerVertexOnFace1 = face1->vertex_index( face1->get_vertex_opposite_to_edge( edge.getID() ) );

      Point3D face1d0 = dx;
      Point3D face1d2 = ( face1->getCoordinates()[outerVertexOnFace1] - x0 ) / ( walberla::real_c( rowsize - 1 ) );

      x = x0 + 1.0 / 3.0 * ( face1d0 + face1d2 ) + dx;
      for ( size_t i = 1; i < rowsize - 2; ++i )
      {
         for ( size_t k = 0; k < srcPtr.size(); ++k )
         {
            srcVector[k] = srcPtr[k][facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )];
         }

         face1->getGeometryMap()->evalF( x, xBlend );
         edgeMemory[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )] =
             expr( xBlend, srcVector );
         x += dx;
      }
   }
}

template < typename ValueType >
inline void assign( const uint_t&                                                              Level,
                    Edge&                                                                      edge,
                    const std::vector< ValueType >&                                            scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                dstId )
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( Level );

   auto dst = edge.getData( dstId )->getPointer( Level );

   // gray south cells
   for ( size_t i = 1; i < rowsize - 2; ++i )
   {
      uint_t    cellIndex = facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE );
      ValueType tmp       = scalars[0] * edge.getData( srcIds[0] )->getPointer( Level )[cellIndex];
      for ( uint_t k = 1; k < srcIds.size(); ++k )
      {
         tmp += scalars[k] * edge.getData( srcIds[k] )->getPointer( Level )[cellIndex];
      }
      dst[cellIndex] = tmp;
   }

   if ( edge.getNumNeighborFaces() == 2 )
   {
      for ( size_t i = 1; i < rowsize - 2; ++i )
      {
         uint_t    cellIndex = facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE );
         ValueType tmp       = scalars[0] * edge.getData( srcIds[0] )->getPointer( Level )[cellIndex];
         for ( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp += scalars[k] * edge.getData( srcIds[k] )->getPointer( Level )[cellIndex];
         }
         dst[cellIndex] = tmp;
      }
   }
}

template < typename ValueType >
inline void multElementwise( const uint_t&                                                              level,
                             Edge&                                                                      edge,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                dstId )
{
   ValueType*                dstPtr = edge.getData( dstId )->getPointer( level );
   std::vector< ValueType* > srcPtr;
   for ( auto src : srcIds )
   {
      srcPtr.push_back( edge.getData( src )->getPointer( level ) );
   }

   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   // gray south cells
   for ( uint_t i = 1; i < rowsize - 2; ++i )
   {
      uint_t    index = facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_SE );
      ValueType tmp   = srcPtr[0][index];
      for ( uint_t k = 1; k < srcIds.size(); ++k )
      {
         tmp *= srcPtr[k][index];
      }
      dstPtr[index] = tmp;
   }

   if ( edge.getNumNeighborFaces() == 2 )
   {
      for ( uint_t i = 1; i < rowsize - 2; ++i )
      {
         uint_t    index = facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_NE );
         ValueType tmp   = srcPtr[0][index];
         for ( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp *= srcPtr[k][index];
         }
         dstPtr[index] = tmp;
      }
   }
}

inline void printFunctionMemory( Edge& edge, const PrimitiveDataID< FunctionMemory< real_t >, Edge >& memoryId, uint_t level )
{
   using namespace std;
   using namespace facedof::macroedge;
   typedef stencilDirection sD;

   uint_t  v_perEdge = hyteg::levelinfo::num_microvertices_per_edge( level );
   real_t* edgeData  = edge.getData( memoryId )->getPointer( level );
   cout << setfill( '=' ) << setw( 100 ) << std::left << "" << endl;
   cout << edge << " South Face ID: " << edge.neighborFaces()[0].getID(); // edge->neighborFaces()[0]->getID().getID();
   if ( edge.getNumHigherDimNeighbors() == 2 )
   {
      cout << " North Face ID: " << edge.neighborFaces()[1].getID();
   }
   cout << setprecision( 6 ) << endl;
   if ( edge.getNumHigherDimNeighbors() == 2 )
   {
      for ( size_t i = 0; i < v_perEdge - 2; ++i )
      {
         cout << setw( 8 ) << setfill( '-' ) << "x"; //edgeData[edge_index(level, i, VERTEX_N)];
      }
      //cout << edgeData[edge_index(level, v_perEdge - 2, VERTEX_N)] << endl << setfill(' ');
      cout << "x" << endl << setfill( ' ' );
      for ( size_t i = 0; i < v_perEdge - 2; ++i )
      {
         cout << "|  \\    ";
      }
      cout << "|  \\" << endl;
      for ( size_t i = 0; i < v_perEdge - 2; ++i )
      {
         cout << "|" << setw( 3 ) << edgeData[indexFaceFromVertex( level, i, sD::CELL_GRAY_NE )];
         cout << "\\" << setw( 3 ) << edgeData[indexFaceFromVertex( level, i + 1, sD::CELL_BLUE_NW )];
      }
      cout << "|" << setw( 3 ) << edgeData[indexFaceFromVertex( level, v_perEdge - 1, sD::CELL_GRAY_NW )] << "\\" << endl;
      for ( size_t i = 0; i < v_perEdge - 2; ++i )
      {
         cout << "|    \\  ";
      }
      cout << "|    \\" << endl;
   }
   //middle vertex
   for ( size_t i = 0; i < v_perEdge - 1; ++i )
   {
      cout << setw( 8 ) << setfill( '-' );
      cout << "x"; //edgeData[edge_index(level, i, VERTEX_C)];
   }
   //cout << edgeData[edge_index(level, v_perEdge - 1, VERTEX_C)] << endl;
   cout << "x" << endl;
   //fill
   cout << "   \\    |";
   for ( size_t i = 0; i < v_perEdge - 2; ++i )
   {
      cout << "  \\    |";
   }
   cout << endl;
   //cell South
   cout << "    \\" << setfill( ' ' ) << setw( 3 ) << edgeData[indexFaceFromVertex( level, 0u, sD::CELL_GRAY_SE )] << "|";
   for ( size_t i = 0; i < v_perEdge - 2; ++i )
   {
      cout << setw( 3 ) << edgeData[indexFaceFromVertex( level, i + 1u, sD::CELL_BLUE_SE )];
      cout << "\\" << setw( 3 ) << edgeData[indexFaceFromVertex( level, i + 1u, sD::CELL_GRAY_SE )] << "|";
   }
   cout << "\n     \\  |";
   for ( size_t i = 0; i < v_perEdge - 2; ++i )
   {
      cout << "    \\  |";
   }

   //vertex South
   cout << "\n        ";
   for ( size_t i = 0; i < v_perEdge - 2; ++i )
   {
      cout << setw( 8 ) << setfill( '-' );
      //cout << edgeData[edge_index(level, i, VERTEX_SE)];
      cout << "x";
   }
   //cout << edgeData[edge_index(level, v_perEdge - 2, VERTEX_SE)] << std::endl;
   cout << "x" << std::endl;
   cout << setfill( '=' ) << setw( 100 ) << "" << setfill( ' ' ) << std::endl;
}

inline void
    printFunctionMemory( Vertex& vertex, const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& memoryId, uint_t level )
{
   real_t* vertexData = vertex.getData( memoryId )->getPointer( level );

   std::cout << std::string( 10, '*' );
   std::cout << " Vertex ID: " << vertex.getID().getID();
   std::cout << " Center: " << vertexData[0];
   std::cout << " Memory ID: " << memoryId;
   std::cout << " " << std::string( 10, '*' ) << std::endl;
   std::cout << "Face ID: |"
             << " Cell " << std::endl;
   for ( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
   {
      std::cout << std::left << std::setw( 9 ) << vertex.neighborFaces()[i].getID() << "|" << vertexData[1 + i] << std::endl;
   }
   std::cout << std::string( 100, '*' ) << std::endl;
}

template < typename ValueType >
inline void add( const uint_t&                                               level,
                 Edge&                                                       edge,
                 const ValueType                                             scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId )
{
   ValueType* dstPtr = edge.getData( dstId )->getPointer( level );

   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   // gray south cells
   for ( uint_t i = 1; i < rowsize - 2; ++i )
   {
      uint_t index = facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_SE );
      dstPtr[index] += scalar;
   }

   if ( edge.getNumNeighborFaces() == 2 )
   {
      for ( uint_t i = 1; i < rowsize - 2; ++i )
      {
         uint_t index = facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_NE );
         dstPtr[index] += scalar;
      }
   }
}

template < typename ValueType >
inline ValueType dot( const uint_t&                                               level,
                      Edge&                                                       edge,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& lhsMemoryId,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsMemoryId )
{
   walberla::math::KahanAccumulator< ValueType > scalarProduct;

   ValueType* lhsPtr = edge.getData( lhsMemoryId )->getPointer( level );
   ValueType* rhsPtr = edge.getData( rhsMemoryId )->getPointer( level );

   size_t rowsize = levelinfo::num_microvertices_per_edge( level );

   // gray south cells
   for ( uint_t i = 1; i < rowsize - 2; ++i )
   {
      uint_t index = facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_SE );
      scalarProduct += lhsPtr[index] * rhsPtr[index];
   }

   if ( edge.getNumNeighborFaces() == 2 )
   {
      // gray north cells
      for ( uint_t i = 1; i < rowsize - 2; ++i )
      {
         uint_t index = facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_NE );
         scalarProduct += lhsPtr[index] * rhsPtr[index];
      }
   }

   return scalarProduct.get();
}

template < typename ValueType >
inline void swap( const uint_t&                                               level,
                  Edge&                                                       edge,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstID )
{
   edge.getData( srcID )->swap( *edge.getData( dstID ), level );
}

} // namespace macroedge
} // namespace facedof
} //namespace hyteg
