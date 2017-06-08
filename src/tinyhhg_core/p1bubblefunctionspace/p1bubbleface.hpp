#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/p1bubblefunctionspace/p1bubblememory.hpp"
#include "p1bubblefacecomm.hpp"

namespace hhg
{
namespace P1BubbleFace
{
//FIXME this can be removed after we are in waberla namespace
using namespace walberla::mpistubs;

namespace CoordsVertex {
enum DirVertex {
  VERTEX_S  = 0,
  VERTEX_SE = 1,
  VERTEX_W  = 2,
  VERTEX_C  = 3,
  VERTEX_E  = 4,
  VERTEX_NW = 5,
  VERTEX_N  = 6,
  CELL_GRAY_SE = 7,
  CELL_GRAY_NW = 8,
  CELL_GRAY_NE = 9,
  CELL_BLUE_SW = 10,
  CELL_BLUE_SE = 11,
  CELL_BLUE_NW = 12
};

const DirVertex neighbors_with_center[13] =
    {VERTEX_C,
     VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W,
     CELL_GRAY_SE, CELL_GRAY_NE, CELL_GRAY_NW,
     CELL_BLUE_SE, CELL_BLUE_NW, CELL_BLUE_SW};
const DirVertex neighbors[12] =
    {VERTEX_S, VERTEX_SE, VERTEX_E, VERTEX_N, VERTEX_NW, VERTEX_W,
     CELL_GRAY_SE, CELL_GRAY_NE, CELL_GRAY_NW,
     CELL_BLUE_SE, CELL_BLUE_NW, CELL_BLUE_SW};

template<size_t Level>
inline size_t index(size_t row, size_t col, DirVertex dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  const size_t grayBaseLength = vertexBaseLength -1;
  const size_t blueBaseLength = vertexBaseLength -2;
  const size_t totalVertices = vertexBaseLength * (vertexBaseLength + 1) / 2;
  const size_t totalCellGray = grayBaseLength * (grayBaseLength + 1) / 2;
  const size_t center = (totalVertices - (vertexBaseLength-row)*(vertexBaseLength-row+1)/2) + col;
  const size_t cellGrayNE = center + totalVertices - row;
  const size_t cellBlueNW = cellGrayNE + (totalCellGray - row) -1;
  switch (dir) {
    case VERTEX_C:
      return center;
    case VERTEX_N:
      return center + vertexBaseLength - row;
    case VERTEX_E:
      return center + 1;
    case VERTEX_S:
      return center - vertexBaseLength - 1 + row;
    case VERTEX_W:
      return center - 1;
    case VERTEX_SE:
      return center - vertexBaseLength + row;
    case VERTEX_NW:
      return center + vertexBaseLength - row - 1;
    case CELL_GRAY_SE:
      return cellGrayNE - (grayBaseLength - row) -1;
    case CELL_GRAY_NE:
      return cellGrayNE;
    case CELL_GRAY_NW:
      return cellGrayNE - 1;
    case CELL_BLUE_SE:
      return cellBlueNW - (blueBaseLength - row);
    case CELL_BLUE_NW:
      return cellBlueNW;
    case CELL_BLUE_SW:
      return cellBlueNW - (blueBaseLength - row) -1;
  }
  return std::numeric_limits<size_t>::max();
}
}//namespace CoordsVertex

namespace CoordsCellGray {
enum Dir {
  VERTEX_SW = 0,
  VERTEX_SE = 1,
  VERTEX_NW = 2,
  CELL_GRAY_C = 3
};

const Dir neighbors[3] = {VERTEX_SE,VERTEX_NW,VERTEX_SW};
const Dir neighbors_with_center[4] = {CELL_GRAY_C,VERTEX_SE,VERTEX_NW,VERTEX_SW};


template<size_t Level>
inline size_t index(size_t row, size_t col, Dir dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  const size_t totalVertices = vertexBaseLength * (vertexBaseLength + 1) / 2;
  const size_t grayBaseLength = vertexBaseLength -1;
  const size_t totalGray = grayBaseLength * (grayBaseLength + 1) / 2;
  const size_t center = totalVertices + totalGray - (grayBaseLength - row) * (grayBaseLength - row  + 1) / 2 + col;
  switch(dir){
    case CELL_GRAY_C:
      return center;
    case VERTEX_SE:
      return center - totalVertices + row + 1;
    case VERTEX_SW:
      return center - totalVertices + row;
    case VERTEX_NW:
      return center - totalVertices + row + vertexBaseLength - row;
  }
  return std::numeric_limits<size_t>::max();
}


}//namesapce CoordsCellGray

namespace CoordsCellBlue {
enum Dir {
  VERTEX_SE = 0,
  VERTEX_NW = 1,
  VERTEX_NE = 2,
  CELL_BLUE_C = 3
};

const Dir neighbors[3] = {VERTEX_SE,VERTEX_NE,VERTEX_NW};
const Dir neighbors_with_center[4] = {CELL_BLUE_C,VERTEX_SE,VERTEX_NE,VERTEX_NW};


template<size_t Level>
inline size_t index(size_t row, size_t col, Dir dir) {
  const size_t vertexBaseLength = levelinfo::num_microvertices_per_edge(Level);
  const size_t totalVertices = vertexBaseLength * (vertexBaseLength + 1) / 2;
  const size_t grayBaseLength = vertexBaseLength -1;
  const size_t totalGray = grayBaseLength * (grayBaseLength + 1) / 2;
  const size_t blueBaseLength = vertexBaseLength -2;
  const size_t totalBlue = blueBaseLength * (blueBaseLength + 1) / 2;
  const size_t center = totalVertices + totalGray + totalBlue - (blueBaseLength - row) * (blueBaseLength - row  + 1) / 2 + col;
  switch(dir){
    case CELL_BLUE_C:
      return center;
    case VERTEX_SE:
      return center - totalVertices + row - totalGray + row + 1;
    case VERTEX_NE:
      return center - totalVertices + row - totalGray + row + vertexBaseLength - row + 1;
    case VERTEX_NW:
      return center - totalVertices + row - totalGray + row + vertexBaseLength - row;
  }
  return std::numeric_limits<size_t>::max();
}

}//namespace CoordsCellBlue


inline void allocate(Face& face, size_t memory_id, size_t minLevel, size_t maxLevel)
{
  face.memory.push_back(new FaceP1Memory());

  for (size_t level = minLevel; level <= maxLevel; ++level)
  {
    getFaceP1BubbleFunctionMemory(face, memory_id)->addlevel(level);
  }
}

inline void free(Face& face, size_t memory_id)
{
  delete face.memory[memory_id];
  face.memory[memory_id] = nullptr;
}

inline void interpolate(Face& face, size_t memory_id, std::function<real_t(const hhg::Point3D&)>& expr, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  Point3D x, x0;

  if (face.edge_orientation[0] == 1)
  {
    x0 = face.edges[0]->v0->coords;
  }
  else
  {
    x0 = face.edges[0]->v1->coords;
  }

  Point3D d0 = face.edge_orientation[0] * face.edges[0]->direction / (walberla::real_c(rowsize-1));
  Point3D d2 = -face.edge_orientation[2] * face.edges[2]->direction / (walberla::real_c(rowsize-1));

  size_t mr_c = 1 + rowsize;
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize-3; ++i)
  {
    x = x0;
    x += (i+1) * d2 + d0;

    for (size_t j = 0; j < inner_rowsize-3; ++j)
    {
      getFaceP1BubbleFunctionMemory(face, memory_id)->data[level][mr_c] = expr(x);
      x += d0;
      mr_c += 1;
    }

    mr_c += 2;
    inner_rowsize -= 1;
  }
}

inline void pull_edges(Face& face, size_t memory_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  real_t* edge_data_0 = NULL;
  real_t* edge_data_1 = NULL;
  real_t* edge_data_2 = NULL;

  MPI_Request req0;
  MPI_Request req1;
  MPI_Request req2;

  int rk = walberla::mpi::MPIManager::instance()->rank();

  if (face.edges[0]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_0 = getEdgeP1Memory(*face.edges[0], memory_id)->data[level];
    }
    else
    {
      MPI_Send(&getEdgeP1Memory(*face.edges[0], memory_id)->data[level][0], rowsize, walberla::MPITrait< real_t >::type(), face.rank, face.edges[0]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_0 = new real_t[rowsize];
    MPI_Irecv(edge_data_0, rowsize, walberla::MPITrait< real_t >::type(), face.edges[0]->rank, face.edges[0]->id, MPI_COMM_WORLD, &req0);
  }

  if (face.edges[1]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_1 = getEdgeP1Memory(*face.edges[1], memory_id)->data[level];
    }
    else
    {
      MPI_Send(&getEdgeP1Memory(*face.edges[1], memory_id)->data[level][0], rowsize, walberla::MPITrait< real_t >::type(), face.rank, face.edges[1]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_1 = new real_t[rowsize];
    MPI_Irecv(edge_data_1, rowsize, walberla::MPITrait< real_t >::type(), face.edges[1]->rank, face.edges[1]->id, MPI_COMM_WORLD, &req1);
  }

  if (face.edges[2]->rank == rk)
  {
    if (face.rank == rk)
    {
      edge_data_2 = getEdgeP1Memory(*face.edges[2], memory_id)->data[level];
    }
    else
    {
      MPI_Send(&getEdgeP1Memory(*face.edges[2], memory_id)->data[level][0], rowsize, walberla::MPITrait< real_t >::type(), face.rank, face.edges[2]->id, MPI_COMM_WORLD);
    }
  }
  else if (face.rank == rk)
  {
    edge_data_2 = new real_t[rowsize];
    MPI_Irecv(edge_data_2, rowsize, walberla::MPITrait< real_t >::type(), face.edges[2]->rank, face.edges[2]->id, MPI_COMM_WORLD, &req2);
  }

  if (face.rank == rk)
  {
    real_t* face_data = getFaceP1BubbleFunctionMemory(face, memory_id)->data[level];

    if (face.edges[0]->rank != rk)
    {
      MPI_Wait(&req0, MPI_STATUS_IGNORE);
    }

    if (face.edges[1]->rank != rk)
    {
      MPI_Wait(&req1, MPI_STATUS_IGNORE);
    }

    if (face.edges[2]->rank != rk)
    {
      MPI_Wait(&req2, MPI_STATUS_IGNORE);
    }

    // edge 0
    if (face.edge_orientation[0] == 1)
    {
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[i] = edge_data_0[i];
      }
    }
    else
    {
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[i] = edge_data_0[rowsize - 1 - i];
      }
    }

    if (face.edges[0]->rank != rk)
    {
      delete[] edge_data_0;
    }

    // edge 1
    if (face.edge_orientation[1] == 1)
    {
      size_t idx = rowsize - 1;
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[idx] = edge_data_1[i];
        idx += rowsize - 1 - i;
      }
    }
    else
    {
      size_t idx = levelinfo::num_microvertices_per_face(level) - 1;
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[idx] = edge_data_1[i];
        idx -= i + 1;
      }
    }

    if (face.edges[1]->rank != rk)
    {
      delete[] edge_data_1;
    }

    // edge 2
    if (face.edge_orientation[2] == 1)
    {
      size_t idx = levelinfo::num_microvertices_per_face(level) - 1;
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[idx] = edge_data_2[i];
        idx -= i+2;
      }
    }
    else
    {
      size_t idx = 0;
      for (size_t i = 0; i < rowsize; ++i)
      {
        face_data[idx] = edge_data_2[i];
        idx += rowsize-i;
      }
    }

    if (face.edges[2]->rank != rk)
    {
      delete[] edge_data_2;
    }
  }
}

inline void assign(Face& face, const std::vector<real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      real_t tmp = scalars[0] * getFaceP1BubbleFunctionMemory(face, src_ids[0])->data[level][mr];

      for (size_t k = 1; k < src_ids.size(); ++k)
      {
        tmp += scalars[k] * getFaceP1BubbleFunctionMemory(face, src_ids[k])->data[level][mr];
      }
      getFaceP1BubbleFunctionMemory(face, dst_id)->data[level][mr] = tmp;

      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }
}

inline void add(Face& face, const std::vector<real_t>& scalars, const std::vector<size_t>& src_ids, size_t dst_id, size_t level)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      real_t tmp = 0.0;

      for (size_t k = 0; k < src_ids.size(); ++k)
      {
        tmp += scalars[k] * getFaceP1BubbleFunctionMemory(face, src_ids[k])->data[level][mr];
      }

      getFaceP1BubbleFunctionMemory(face, dst_id)->data[level][mr] += tmp;

      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }
}

inline real_t dot(Face& face, size_t lhs_id, size_t rhs_id, size_t level)
{
  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  size_t inner_rowsize = rowsize;

  size_t mr = 1 + rowsize;

  for (size_t i = 0; i < rowsize - 3; ++i)
  {
    for (size_t j = 0; j < inner_rowsize - 3; ++j)
    {
      sp += getFaceP1BubbleFunctionMemory(face, lhs_id)->data[level][mr] * getFaceP1BubbleFunctionMemory(face, rhs_id)->data[level][mr];
      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }

  return sp;
}

template<size_t Level>
inline void apply_tmpl(Face& face, size_t opr_id, size_t src_id, size_t dst_id, UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto opr_data = getFaceP1BubbleStencilMemory(face, opr_id)->data[Level];

  real_t* face_vertex_stencil = opr_data[0];
  real_t* face_gray_stencil = opr_data[1];
  real_t* face_blue_stencil = opr_data[2];

  real_t* src = getFaceP1BubbleFunctionMemory(face, src_id)->data[Level];
  real_t* dst = getFaceP1BubbleFunctionMemory(face, dst_id)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize - 2; ++i)
  {
    for (size_t j = 1; j  < inner_rowsize - 2; ++j)
    {
      tmp = face_vertex_stencil[CoordsVertex::VERTEX_C] * src[CoordsVertex::index<Level>(i, j, CoordsVertex::VERTEX_C)];

      for (auto neighbor : CoordsVertex::neighbors)
      {
        tmp += face_vertex_stencil[neighbor] * src[CoordsVertex::index<Level>(i, j, neighbor)];
      }

      if (update == Replace) {
        dst[CoordsVertex::index<Level>(i, j, CoordsVertex::VERTEX_C)] = tmp;
      } else if (update == Add) {
        dst[CoordsVertex::index<Level>(i, j, CoordsVertex::VERTEX_C)] += tmp;
      }
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      // TODO: how to do this better?
      if ((i == 0 && j == 0) || (i == 0 && j == rowsize - 2) || (i == rowsize - 2 && j == 0)) {
        continue;
      }

      tmp = face_gray_stencil[CoordsCellGray::CELL_GRAY_C] * src[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)];

      for (auto neighbor : CoordsCellGray::neighbors)
      {
        tmp += face_gray_stencil[neighbor] * src[CoordsCellGray::index<Level>(i, j, neighbor)];
      }

      if (update == Replace) {
        dst[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)] = tmp;
      } else if (update == Add) {
        dst[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)] += tmp;
      }
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      tmp = face_blue_stencil[CoordsCellBlue::CELL_BLUE_C] * src[CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)];

      for (auto neighbor : CoordsCellBlue::neighbors)
      {
        tmp += face_blue_stencil[neighbor] * src[CoordsCellBlue::index<Level>(i, j, neighbor)];
      }

      if (update == Replace) {
        dst[CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)] = tmp;
      } else if (update == Add) {
        dst[CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, apply_tmpl, apply)

//template<size_t Level>
//inline void smooth_gs_tmpl(Face& face, size_t opr_id, size_t dst_id, size_t rhs_id)
//{
//  using namespace CoordsVertex;
//  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
//  size_t inner_rowsize = rowsize;
//
//  real_t* opr_data = getFaceStencilMemory(face, opr_id)->data[Level];
//  real_t* dst = getFaceP1BubbleFunctionMemory(face, dst_id)->data[Level];
//  real_t* rhs = getFaceP1BubbleFunctionMemory(face, rhs_id)->data[Level];
//
//  real_t tmp;
//
//  for (size_t i = 1; i < rowsize - 2; ++i)
//  {
//    for (size_t j = 1; j  < inner_rowsize - 2; ++j)
//    {
//      tmp = rhs[index<Level>(i, j, VERTEX_C)];
//
//      for (auto neighbor : neighbors)
//      {
//        tmp -= opr_data[neighbor] * dst[index<Level>(i, j, neighbor)];
//      }
//
//      dst[index<Level>(i, j, VERTEX_C)] = tmp / opr_data[VERTEX_C];
//    }
//    --inner_rowsize;
//  }
//}
//
//SPECIALIZE(void, smooth_gs_tmpl, smooth_gs)
//
//template<size_t Level>
//inline void prolongate_tmpl(Face& face, size_t memory_id)
//{
//  using namespace CoordsVertex;
//  size_t N_c = levelinfo::num_microvertices_per_edge(Level);
//  size_t N_c_i = N_c;
//
//  real_t* v_f = getFaceP1BubbleFunctionMemory(face, memory_id)->data[Level+1];
//  real_t* v_c = getFaceP1BubbleFunctionMemory(face, memory_id)->data[Level];
//
//  size_t j;
//
//  for (size_t i = 1; i < N_c-1; ++i)
//  {
//    for (j = 1; j < N_c_i-2; ++j)
//    {
//      v_f[index<Level+1>(2*i, 2*j, VERTEX_C)] = v_c[index<Level>(i, j, VERTEX_C)];
//      v_f[index<Level+1>(2*i - 1, 2*j - 1, VERTEX_C)] = 0.5 * (v_c[index<Level>(i-1, j, VERTEX_C)] + v_c[index<Level>(i, j-1, VERTEX_C)]);
//      v_f[index<Level+1>(2*i - 1, 2*j, VERTEX_C)] = 0.5 * (v_c[index<Level>(i, j, VERTEX_C)] + v_c[index<Level>(i-1, j, VERTEX_C)]);
//      v_f[index<Level+1>(2*i, 2*j - 1, VERTEX_C)] = 0.5 * (v_c[index<Level>(i, j, VERTEX_C)] + v_c[index<Level>(i, j-1, VERTEX_C)]);
//    }
//
//    v_f[index<Level+1>(2*i - 1, 2*j - 1, VERTEX_C)] = 0.5 * (v_c[index<Level>(i-1, j, VERTEX_C)] + v_c[index<Level>(i, j-1, VERTEX_C)]);
//    v_f[index<Level+1>(2*i - 1, 2*j, VERTEX_C)] = 0.5 * (v_c[index<Level>(i, j, VERTEX_C)] + v_c[index<Level>(i-1, j, VERTEX_C)]);
//    v_f[index<Level+1>(2*i, 2*j - 1, VERTEX_C)] = 0.5 * (v_c[index<Level>(i, j, VERTEX_C)] + v_c[index<Level>(i, j-1, VERTEX_C)]);
//
//    --N_c_i;
//  }
//}
//
//SPECIALIZE(void, prolongate_tmpl, prolongate)
//
//template<size_t Level>
//inline void restrict_tmpl(Face& face, size_t memory_id)
//{
//  using namespace CoordsVertex;
//  size_t N_c = levelinfo::num_microvertices_per_edge(Level-1);
//  size_t N_c_i = N_c;
//
//  real_t* v_f = getFaceP1BubbleFunctionMemory(face, memory_id)->data[Level];
//  real_t* v_c = getFaceP1BubbleFunctionMemory(face, memory_id)->data[Level-1];
//
//  real_t tmp;
//
//  for (size_t i = 1; i < N_c - 2; ++i)
//  {
//    for (size_t j = 1; j < N_c_i - 2; ++j)
//    {
//      tmp = v_f[index<Level>(2*i, 2*j, VERTEX_C)];
//
//      for (auto neighbor : neighbors)
//      {
//        tmp += 0.5 * v_f[index<Level>(2*i, 2*j, neighbor)];
//      }
//
//      v_c[index<Level-1>(i, j, VERTEX_C)] = tmp;
//    }
//
//    --N_c_i;
//  }
//}
//
//SPECIALIZE(void, restrict_tmpl, restrict)
//
///// Checks if a given index is a the boundary of the face
///// \param index The index which should be checked
///// \param length Size of the triangle in the first dimension
//bool is_boundary(size_t index, size_t length)
//{
//  if(index < length) return true;
//  while(index >= length){
//    index -= length;
//    length--;
//  }
//  return(index == 0 || index == (length -1));
//}
//
//inline void printmatrix(Face& face, size_t opr_id, size_t src_id, size_t level)
//{
//  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
//  size_t inner_rowsize = rowsize;
//
//  real_t* opr_data = getFaceStencilMemory(face, opr_id)->data[level];
//  real_t* src = getFaceP1BubbleFunctionMemory(face, src_id)->data[level];
//  size_t br = 1;
//  size_t mr = 1 + rowsize ;
//  size_t tr = mr + (rowsize - 1);
//
//  for (size_t i = 0; i < rowsize - 3; ++i)
//  {
//    for (size_t j = 0; j < inner_rowsize - 3; ++j)
//    {
//      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[br], opr_data[0]);
//      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[br+1], opr_data[1]);
//      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[mr-1], opr_data[2]);
//      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[mr], opr_data[3]);
//      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[mr+1], opr_data[4]);
//      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[tr-1], opr_data[5]);
//      fmt::printf("%d\t%d\t%e\n", (size_t)src[mr], (size_t)src[tr], opr_data[6]);
//
//      br += 1;
//      mr += 1;
//      tr += 1;
//    }
//
//    br += 3;
//    mr += 2;
//    tr += 1;
//    --inner_rowsize;
//  }
//}

}// namespace P1BubbleFace
}// namespace hhg
