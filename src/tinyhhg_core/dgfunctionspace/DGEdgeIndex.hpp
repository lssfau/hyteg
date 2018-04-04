#pragma once

#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleEdgeIndex.hpp"

namespace hhg {
namespace DGEdge {

using walberla::uint_t;

constexpr inline uint_t indexDGFaceFromVertex( const uint_t& level, const uint_t pos, stencilDirection dir )
{
   return hhg::BubbleEdge::indexFaceFromVertex( level, pos, dir );
}

inline void printFunctionMemory( Edge& edge, const PrimitiveDataID< FunctionMemory< real_t >, Edge >& memoryId, uint_t level )
{
   using namespace std;
   using namespace hhg::BubbleEdge;
   typedef stencilDirection sD;

   uint_t  v_perEdge = hhg::levelinfo::num_microvertices_per_edge( level );
   real_t* edgeData  = edge.getData( memoryId )->getPointer( level );
   cout << setfill( '=' ) << setw( 100 ) << std::left << "" << endl;
   cout << edge << " South Face ID: " << edge.neighborFaces()[0].getID(); // edge->neighborFaces()[0]->getID().getID();
   if( edge.getNumHigherDimNeighbors() == 2 )
   {
      cout << " North Face ID: " << edge.neighborFaces()[1].getID();
   }
   cout << setprecision( 6 ) << endl;
   if( edge.getNumHigherDimNeighbors() == 2 )
   {
      for( size_t i = 0; i < v_perEdge - 2; ++i )
      {
         cout << setw( 8 ) << setfill( '-' ) << "x"; //edgeData[edge_index(level, i, VERTEX_N)];
      }
      //cout << edgeData[edge_index(level, v_perEdge - 2, VERTEX_N)] << endl << setfill(' ');
      cout << "x" << endl << setfill( ' ' );
      for( size_t i = 0; i < v_perEdge - 2; ++i )
      {
         cout << "|  \\    ";
      }
      cout << "|  \\" << endl;
      for( size_t i = 0; i < v_perEdge - 2; ++i )
      {
         cout << "|" << setw( 3 ) << edgeData[indexFaceFromVertex( level, i, sD::CELL_GRAY_NE )];
         cout << "\\" << setw( 3 ) << edgeData[indexFaceFromVertex( level, i + 1, sD::CELL_BLUE_NW )];
      }
      cout << "|" << setw( 3 ) << edgeData[indexFaceFromVertex( level, v_perEdge - 1, sD::CELL_GRAY_NW )] << "\\" << endl;
      for( size_t i = 0; i < v_perEdge - 2; ++i )
      {
         cout << "|    \\  ";
      }
      cout << "|    \\" << endl;
   }
   //middle vertex
   for( size_t i = 0; i < v_perEdge - 1; ++i )
   {
      cout << setw( 8 ) << setfill( '-' );
      cout << "x"; //edgeData[edge_index(level, i, VERTEX_C)];
   }
   //cout << edgeData[edge_index(level, v_perEdge - 1, VERTEX_C)] << endl;
   cout << "x" << endl;
   //fill
   cout << "   \\    |";
   for( size_t i = 0; i < v_perEdge - 2; ++i )
   {
      cout << "  \\    |";
   }
   cout << endl;
   //cell South
   cout << "    \\" << setfill( ' ' ) << setw( 3 ) << edgeData[indexFaceFromVertex( level, 0u, sD::CELL_GRAY_SE )] << "|";
   for( size_t i = 0; i < v_perEdge - 2; ++i )
   {
      cout << setw( 3 ) << edgeData[indexFaceFromVertex( level, i + 1u, sD::CELL_BLUE_SE )];
      cout << "\\" << setw( 3 ) << edgeData[indexFaceFromVertex( level, i + 1u, sD::CELL_GRAY_SE )] << "|";
   }
   cout << "\n     \\  |";
   for( size_t i = 0; i < v_perEdge - 2; ++i )
   {
      cout << "    \\  |";
   }

   //vertex South
   cout << "\n        ";
   for( size_t i = 0; i < v_perEdge - 2; ++i )
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
   for( uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i )
   {
      std::cout << std::left << std::setw( 9 ) << vertex.neighborFaces()[i].getID() << "|" << vertexData[1 + i] << std::endl;
   }
   std::cout << std::string( 100, '*' ) << std::endl;
}

} // namespace DGEdge
} // namespace hhg