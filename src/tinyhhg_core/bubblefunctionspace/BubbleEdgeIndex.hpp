#pragma once

#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleMemory.hpp"

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

constexpr std::array<DirVertex,6> neighbors =
  {{CELL_GRAY_SE, CELL_GRAY_NE, CELL_GRAY_NW, CELL_GRAY_SW,
    CELL_BLUE_SE, CELL_BLUE_NW}};

constexpr std::array<DirVertex,3> neighbors_south =
  {{CELL_GRAY_SW, CELL_BLUE_SE, CELL_GRAY_SE}};

constexpr std::array<DirVertex,3> neighbors_north =
  {{CELL_GRAY_NW, CELL_BLUE_NW, CELL_GRAY_NE}};

//first face is south face by convention

template<size_t Level>
constexpr inline size_t index(size_t pos, DirVertex dir) {
  const size_t vertexOnEdge = levelinfo::num_microvertices_per_edge(Level);
  //WALBERLA_ASSERT_LESS_EQUAL(pos,vertexOnEdge);
  const size_t startFaceS = 0;
  const size_t startFaceN = 2 * (vertexOnEdge - 1) - 1;
  switch (dir) {
    case CELL_GRAY_SE:
      return startFaceS + pos * 2;
    case CELL_GRAY_NE:
      return startFaceN + pos * 2;
    case CELL_GRAY_NW:
      return startFaceN + pos * 2 - 2;
    case CELL_GRAY_SW:
      return startFaceS + (pos -1) * 2;
    case CELL_BLUE_SE:
      return startFaceS + pos * 2 -1;
    case CELL_BLUE_NW:
      return startFaceN + pos * 2 - 1;
    default:
      //WALBERLA_ASSERT(false, "wrong dir");
      return std::numeric_limits<size_t>::max();
  }

}

SPECIALIZE(size_t, index, edge_index)

}//namespace EdgeCoordsVertex

inline void printFunctionMemory(Edge &edge, const PrimitiveDataID<EdgeBubbleFunctionMemory< real_t >, Edge> &memoryId, uint_t level)
{
  using namespace std;
  using namespace hhg::BubbleEdge::EdgeCoordsVertex;

  uint_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(level);
  real_t* edgeData = edge.getData(memoryId)->data[level].get();
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
           << edgeData[edge_index(level, i, CELL_GRAY_NE)];
      cout << "\\" << setw(3) << edgeData[edge_index(level, i + 1, CELL_BLUE_NW)];
    }
    cout << "|" << setw(3) << edgeData[edge_index(level, v_perEdge - 1, CELL_GRAY_NW)] << "\\" << endl;
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
  cout << "    \\" << setfill(' ') << setw(3) << edgeData[edge_index(level, 0u, CELL_GRAY_SE)] << "|";
  for (size_t i = 0; i < v_perEdge - 2; ++i)
  {
    cout << setw(3) << edgeData[edge_index(level, i + 1u, CELL_BLUE_SE)];
    cout << "\\" << setw(3) << edgeData[edge_index(level, i + 1u, CELL_GRAY_SE)] << "|";
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
  real_t* vertexData = vertex.getData(memoryId)->data[level].get();

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
