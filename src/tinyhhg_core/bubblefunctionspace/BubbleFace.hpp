#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "BubbleMemory.hpp"
#include "BubbleFaceIndex.hpp"

namespace hhg {
namespace BubbleFace {

template<size_t Level>
inline void assign_tmpl(Face &face,
                        const std::vector<real_t> &scalars,
                        const std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face>> &srcIds,
                        const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &dstId) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      real_t tmp = scalars[0] * face.getData(srcIds[0])->data[Level][CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)];

      for (size_t k = 1; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->data[Level][CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)];
      }
      face.getData(dstId)->data[Level][CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)] = tmp;
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      real_t tmp = scalars[0] * face.getData(srcIds[0])->data[Level][CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)];

      for (size_t k = 1; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->data[Level][CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)];
      }
      face.getData(dstId)->data[Level][CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, assign_tmpl, assign)

template<size_t Level>
inline void add_tmpl(Face &face,
                     const std::vector<real_t> &scalars,
                     const std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face>> &srcIds,
                     const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &dstId) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      real_t tmp = 0.0;

      for (size_t k = 0; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->data[Level][CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)];
      }
      face.getData(dstId)->data[Level][CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)] += tmp;
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      real_t tmp = 0.0;

      for (size_t k = 0; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->data[Level][CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)];
      }
      face.getData(dstId)->data[Level][CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)] += tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, add_tmpl, add)

template<size_t Level>
inline real_t dot_tmpl(Face &face,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &lhsId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &rhsId) {
  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto &lhs_data = face.getData(lhsId)->data[Level];
  auto &rhs_data = face.getData(rhsId)->data[Level];

  for (size_t i = 0; i < rowsize - 1; ++i) {
    for (size_t j = 0; j < inner_rowsize - 1; ++j) {
      sp += lhs_data[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)]
          *rhs_data[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)];
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i) {
    for (size_t j = 0; j < inner_rowsize - 2; ++j) {
      sp += lhs_data[CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)]
          *rhs_data[CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)];
    }
    --inner_rowsize;
  }

  return sp;
}

SPECIALIZE(real_t, dot_tmpl, dot)

template<size_t Level>
inline void apply_tmpl(Face& face, const PrimitiveDataID<FaceBubbleStencilMemory, Face>& operatorId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &srcId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &dstId, UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = face.getData(operatorId)->data[Level];

  auto& face_gray_stencil = opr_data[0];
  auto& face_blue_stencil = opr_data[1];

  auto& src = face.getData(srcId)->data[Level];
  auto& dst = face.getData(dstId)->data[Level];

  real_t tmp;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      tmp = face_gray_stencil[CoordsCellGray::CELL_GRAY_C] * src[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)];

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

template<size_t Level>
inline void enumerate_tmpl(Face &face, const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &dstId, uint_t& num) {
  using walberla::real_c;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      face.getData(dstId)->data[Level][CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)] = real_c(num++);
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      face.getData(dstId)->data[Level][CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)] = real_c(num++);
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, enumerate_tmpl, enumerate)

template<size_t Level>
inline void saveOperator_tmpl(Face& face, const PrimitiveDataID<FaceBubbleStencilMemory, Face>& operatorId,
                              const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &srcId,
                              const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &dstId, std::ostream& out)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = face.getData(operatorId)->data[Level];

  auto& face_gray_stencil = opr_data[0];
  auto& face_blue_stencil = opr_data[1];

  auto& src = face.getData(srcId)->data[Level];
  auto& dst = face.getData(dstId)->data[Level];

  real_t tmp;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      tmp = face_gray_stencil[CoordsCellGray::CELL_GRAY_C] * src[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)];

      out << fmt::format("{}\t{}\t{}\n", dst[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)], src[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)], face_gray_stencil[CoordsCellGray::CELL_GRAY_C]);
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      tmp = face_blue_stencil[CoordsCellBlue::CELL_BLUE_C] * src[CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)];

      out << fmt::format("{}\t{}\t{}\n", dst[CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)], src[CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)], face_blue_stencil[CoordsCellBlue::CELL_BLUE_C]);
    }
    --inner_rowsize;
  }
}

SPECIALIZE(void, saveOperator_tmpl, saveOperator)

template<size_t Level>
inline void printFunctionMemory(Face& face, const PrimitiveDataID<FaceBubbleFunctionMemory< real_t >, Face> &dstId){
  using namespace std;
  real_t* faceMemory = face.getData(dstId)->data[Level].get();
  uint_t verticesPerDge = hhg::levelinfo::num_microvertices_per_edge(Level);
  cout << setfill('=') << setw(100) << "" << endl;
  cout << face << std::left << setprecision(1) << fixed << setfill(' ') << endl;
  std::cout << "Cell Blue: " << std::endl;
  for (size_t i = 0; i < verticesPerDge-2; ++i) {
    for (size_t j = 0; j < verticesPerDge-2 - i; ++j) {
      cout << setw(5) << faceMemory[CoordsCellBlue::index<Level>(i, j, CoordsCellBlue::CELL_BLUE_C)] << "|";
    }
    std::cout << std::endl;
  }
  cout << "Cell Gray: " << std::endl;
  for (size_t i = 0; i < verticesPerDge-1; ++i) {
    for (size_t j = 0; j < verticesPerDge-1 - i; ++j) {
      cout << setw(5) << faceMemory[CoordsCellGray::index<Level>(i, j, CoordsCellGray::CELL_GRAY_C)] << "|";
    }
    std::cout << std::endl;
  }
  cout << setw(100) << setfill(' ') << endl;

}

}// namespace P1BubbleFace
}// namespace hhg
