#ifndef P2VERTEX_HPP
#define P2VERTEX_HPP

#include <levelinfo.hpp>
#include <comm.hpp>

#include <fmt/format.h>

namespace hhg
{

/// P2Vertex namespace for P2 macro-vertex kernels
/// 
/// [Vertex dof, Edges: [Edge ghost dof, Vertex ghost dof], Faces: [Edge ghost dof]]
namespace P2Vertex
{

inline void allocate(Vertex& vertex, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  vertex.data.push_back(std::vector<double*>());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    size_t num_deps = 3 * vertex.edges.size();
    size_t total_n_dofs = levelinfo::num_microvertices_per_vertex(level) + num_deps;
    double* new_data = new double[total_n_dofs];
    memset(new_data, 0, total_n_dofs * sizeof(double));
    vertex.data[memory_id].push_back(new_data);
  }
}

inline void free(Vertex& vertex, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    delete[] vertex.data[memory_id][level - minLevel];
  }
}

inline void interpolate(Vertex& vertex, size_t memory_id, std::function<double(const hhg::Point3D&)>& expr, size_t level)
{
  vertex.data[memory_id][level-2][0] = expr(vertex.coords);
}

inline void print(Vertex & vertex, size_t memory_id, size_t level) {
    //bad hack!
//    for(size_t i = 0; i < 10000;++i){
//      fmt::print(" {} ",vertex.data[memory_id][level - 2][i]);
//    }

  auto vertexData = vertex.data[memory_id][level -2];
  int pos = 0;
  fmt::print("Center dof: {}\n" , vertexData[pos++]);
  fmt::print("{:<8} {:<8} {:<8} {:<8}  -Edges-\n","inner","outer","localID","globalID");
  for(size_t i = 0; i < vertex.edges.size(); ++i){
    fmt::print("{:<8} {:<8} {:<8} {:<8}\n",vertexData[pos],vertexData[pos+1],i,vertex.edges[i]->id);
    pos += 2;
  }
  fmt::print("{:<8} {:<8} {:<8}  -Faces-\n","edgeDoF","localID","globalID");
  for(size_t i = 0; i < vertex.faces.size(); ++i){
    fmt::print("{:<8} {:<8} {:<8}\n",vertexData[pos++],i,vertex.faces[i]->id);
    //fmt::print("{}\t{}\t{}\n",vertexData[pos++],i,vertex.faces[i]->id);
  }

}

}
}

#endif /* P2VERTEX_HPP */