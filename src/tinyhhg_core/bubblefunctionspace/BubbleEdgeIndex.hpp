#pragma once

#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleMemory.hpp"
#include "tinyhhg_core/StencilDirections.hpp"

namespace hhg{
namespace BubbleEdge{
//FIXME this can be removed after we moved into walberla namespace
using namespace walberla::mpistubs;

namespace EdgeCoordsVertex {
enum DirVertex {
  CELL_GRAY_SW = 0,
  CELL_BLUE_SE = 1,
  CELL_GRAY_SE = 2,
  CELL_GRAY_NW = 3,
  CELL_BLUE_NW = 4,
  CELL_GRAY_NE = 5
};

}//namespace EdgeCoordsVertex

constexpr inline uint_t indexEdgeStencil(const stencilDirection dir){
  typedef hhg::stencilDirection sD;
  switch(dir) {
    case sD::CELL_GRAY_SW:
      return 0;
    case sD::CELL_BLUE_SE:
      return 1;
    case sD::CELL_GRAY_SE:
      return 2;
    case sD::CELL_GRAY_NW:
      return 3;
    case sD::CELL_BLUE_NW:
      return 4;
    case sD::CELL_GRAY_NE:
      return 5;
    default:
      return std::numeric_limits<size_t>::max();
  }
}

constexpr std::array<stencilDirection,6> neighbors =
  {{stencilDirection::CELL_GRAY_SE, stencilDirection::CELL_GRAY_NE, stencilDirection::CELL_GRAY_NW, stencilDirection::CELL_GRAY_SW,
    stencilDirection::CELL_BLUE_SE, stencilDirection::CELL_BLUE_NW}};

constexpr std::array<stencilDirection,3> neighbors_south =
  {{stencilDirection::CELL_GRAY_SW, stencilDirection::CELL_BLUE_SE, stencilDirection::CELL_GRAY_SE}};

constexpr std::array<stencilDirection,3> neighbors_north =
  {{stencilDirection::CELL_GRAY_NW, stencilDirection::CELL_BLUE_NW, stencilDirection::CELL_GRAY_NE}};

//first face is south face by convention

template<size_t Level>
constexpr inline size_t indexFaceFromVertex(size_t pos, stencilDirection dir) {
  typedef stencilDirection sD;
  const size_t vertexOnEdge = levelinfo::num_microvertices_per_edge(Level);
  WALBERLA_ASSERT_LESS_EQUAL(pos,vertexOnEdge);
  const size_t startFaceS = 0;
  const size_t startFaceN = 2 * (vertexOnEdge - 1) - 1;
  switch (dir) {
    case sD::CELL_GRAY_SE:
      return startFaceS + pos * 2;
    case sD::CELL_GRAY_NE:
      return startFaceN + pos * 2;
    case sD::CELL_GRAY_NW:
      return startFaceN + pos * 2 - 2;
    case sD::CELL_GRAY_SW:
      return startFaceS + (pos -1) * 2;
    case sD::CELL_BLUE_SE:
      return startFaceS + pos * 2 -1;
    case sD::CELL_BLUE_NW:
      return startFaceN + pos * 2 - 1;
    default:
      WALBERLA_ASSERT(false);
      return std::numeric_limits<size_t>::max();
  }

}

SPECIALIZE(size_t, indexFaceFromVertex, edge_index)

inline void printFunctionMemory(Edge &edge, const PrimitiveDataID<EdgeBubbleFunctionMemory< real_t >, Edge> &memoryId, uint_t level)
{
  using namespace std;
  using namespace hhg::BubbleEdge;
  typedef stencilDirection sD;

  uint_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(level);
  real_t* edgeData = edge.getData(memoryId)->getPointer( level );
  cout << setfill('=') << setw(100) << std::left << "" << endl;
  cout << edge << " South Face ID: " << edge.neighborFaces()[0].getID();  // edge->neighborFaces()[0]->getID().getID();
  if (edge.getNumHigherDimNeighbors() == 2)
  { cout << " North Face ID: " << edge.neighborFaces()[1].getID(); }
  cout << setprecision(6) << endl;
  if (edge.getNumHigherDimNeighbors() == 2)
  {
    for (size_t i = 0; i < v_perEdge - 2; ++i)
    {
      cout << setw(8) << setfill('-') << "x"; //edgeData[edge_index(level, i, VERTEX_N)];
    }
    //cout << edgeData[edge_index(level, v_perEdge - 2, VERTEX_N)] << endl << setfill(' ');
    cout << "x" << endl << setfill(' ');
    for (size_t i = 0; i < v_perEdge - 2; ++i)
    { cout << "|  \\    "; }
    cout << "|  \\" << endl;
    for (size_t i = 0; i < v_perEdge - 2; ++i)
    {
      cout << "|" << setw(3)
           << edgeData[edge_index(level, i, sD::CELL_GRAY_NE)];
      cout << "\\" << setw(3) << edgeData[edge_index(level, i + 1, sD::CELL_BLUE_NW)];
    }
    cout << "|" << setw(3) << edgeData[edge_index(level, v_perEdge - 1, sD::CELL_GRAY_NW)] << "\\" << endl;
    for (size_t i = 0; i < v_perEdge - 2; ++i)
    { cout << "|    \\  "; }
    cout << "|    \\" << endl;
  }
//middle vertex
  for (size_t i = 0; i < v_perEdge - 1; ++i)
  {
    cout << setw(8) << setfill('-');
    cout << "x"; //edgeData[edge_index(level, i, VERTEX_C)];
  }
  //cout << edgeData[edge_index(level, v_perEdge - 1, VERTEX_C)] << endl;
  cout << "x" << endl;
//fill
  cout << "   \\    |";
  for (size_t i = 0; i < v_perEdge - 2; ++i)
  { cout << "  \\    |"; }
  cout << endl;
//cell South
  cout << "    \\" << setfill(' ') << setw(3) << edgeData[edge_index(level, 0u, sD::CELL_GRAY_SE)] << "|";
  for (size_t i = 0; i < v_perEdge - 2; ++i)
  {
    cout << setw(3) << edgeData[edge_index(level, i + 1u, sD::CELL_BLUE_SE)];
    cout << "\\" << setw(3) << edgeData[edge_index(level, i + 1u, sD::CELL_GRAY_SE)] << "|";
  }
  cout << "\n     \\  |";
  for (size_t i = 0; i < v_perEdge - 2; ++i)
  { cout << "    \\  |"; }

//vertex South
  cout << "\n        ";
  for (size_t i = 0; i < v_perEdge - 2; ++i)
  {
    cout << setw(8) << setfill('-');
    //cout << edgeData[edge_index(level, i, VERTEX_SE)];
    cout << "x";
  }
  //cout << edgeData[edge_index(level, v_perEdge - 2, VERTEX_SE)] << std::endl;
  cout << "x" << std::endl;
  cout << setfill('=') << setw(100) << "" << setfill(' ') << std::endl;
}

inline void printFunctionMemory(Vertex& vertex, const PrimitiveDataID<VertexBubbleFunctionMemory< real_t >, Vertex> &memoryId, uint_t level)
{
  real_t* vertexData = vertex.getData(memoryId)->getPointer( level );

  std::cout <<  std::string(10,'*');
  std::cout << " Vertex ID: " << vertex.getID().getID();
  std::cout << " Center: " << vertexData[0];
  std::cout << " Memory ID: "<< memoryId;
  std::cout <<  " " << std::string(10,'*') << std::endl;
  std::cout << "Face ID: |" << " Cell " << std::endl;
  for(uint_t i = 0; i < vertex.getNumNeighborFaces(); ++i){
    std::cout << std::left << std::setw(9) << vertex.neighborFaces()[i].getID() << "|"  << vertexData[1+i] << std::endl;
  }
  std::cout <<  std::string(100,'*') << std::endl;
}

}//namespace BubbleToP1Edge
}//namespace hhg
