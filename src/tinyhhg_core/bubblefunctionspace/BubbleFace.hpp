#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "BubbleMemory.hpp"
#include "BubbleFaceIndex.hpp"

namespace hhg {
namespace BubbleFace {

template< typename ValueType, size_t Level >
inline void assign_tmpl(Face &face,
                        const std::vector<ValueType> &scalars,
                        const std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face>> &srcIds,
                        const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      ValueType tmp = scalars[0] * face.getData(srcIds[0])->getPointer( Level )[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)];

      for (size_t k = 1; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->getPointer( Level )[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)];
      }
      face.getData(dstId)->getPointer( Level )[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)] = tmp;
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      ValueType tmp = scalars[0] * face.getData(srcIds[0])->getPointer( Level )[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)];

      for (size_t k = 1; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->getPointer( Level )[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)];
      }
      face.getData(dstId)->getPointer( Level )[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)] = tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, assign_tmpl, assign)

template< typename ValueType, size_t Level >
inline void add_tmpl(Face &face,
                     const std::vector<ValueType> &scalars,
                     const std::vector<PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face>> &srcIds,
                     const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      ValueType tmp = 0.0;

      for (size_t k = 0; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->getPointer( Level )[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)];
      }
      face.getData(dstId)->getPointer( Level )[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)] += tmp;
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      ValueType tmp = 0.0;

      for (size_t k = 0; k < srcIds.size(); ++k)
      {
        tmp += scalars[k] * face.getData(srcIds[k])->getPointer( Level )[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)];
      }
      face.getData(dstId)->getPointer( Level )[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)] += tmp;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, add_tmpl, add)

template< typename ValueType, size_t Level >
inline real_t dot_tmpl(Face &face,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &lhsId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &rhsId) {
  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto lhs_data = face.getData(lhsId)->getPointer( Level );
  auto rhs_data = face.getData(rhsId)->getPointer( Level );

  for (size_t i = 0; i < rowsize - 1; ++i) {
    for (size_t j = 0; j < inner_rowsize - 1; ++j) {
      sp += lhs_data[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)]
          *rhs_data[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)];
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i) {
    for (size_t j = 0; j < inner_rowsize - 2; ++j) {
      sp += lhs_data[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)]
          *rhs_data[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)];
    }
    --inner_rowsize;
  }

  return sp;
}

SPECIALIZE_WITH_VALUETYPE(real_t, dot_tmpl, dot)

template< typename ValueType, size_t Level >
inline void apply_tmpl(Face& face, const PrimitiveDataID<FaceBubbleStencilMemory, Face>& operatorId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &srcId,
                       const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId, UpdateType update)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = face.getData(operatorId)->data[Level];

  auto& face_gray_stencil = opr_data[0];
  auto& face_blue_stencil = opr_data[1];

  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  ValueType tmp;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      tmp = face_gray_stencil[FaceCoordsCellGray::CELL_GRAY_C] * src[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)];

      if (update == Replace) {
        dst[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)] = tmp;
      } else if (update == Add) {
        dst[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)] += tmp;
      }
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      tmp = face_blue_stencil[FaceCoordsCellBlue::CELL_BLUE_C] * src[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)];

      if (update == Replace) {
        dst[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)] = tmp;
      } else if (update == Add) {
        dst[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, apply_tmpl, apply)

template< typename ValueType, size_t Level >
inline void enumerate_tmpl(Face &face, const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId, uint_t& num) {
  using walberla::real_c;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      face.getData(dstId)->getPointer( Level )[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)] = real_c(num++);
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      face.getData(dstId)->getPointer( Level )[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)] = real_c(num++);
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, enumerate_tmpl, enumerate)

template< typename ValueType, size_t Level >
inline void saveOperator_tmpl(Face& face, const PrimitiveDataID<FaceBubbleStencilMemory, Face>& operatorId,
                              const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &srcId,
                              const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId, std::ostream& out)
{
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto& opr_data = face.getData(operatorId)->data[Level];

  auto& face_gray_stencil = opr_data[0];
  auto& face_blue_stencil = opr_data[1];

  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  ValueType tmp;

  for (size_t i = 0; i < rowsize - 1; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 1; ++j)
    {
      tmp = face_gray_stencil[FaceCoordsCellGray::CELL_GRAY_C] * src[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)];

      out << fmt::format("{}\t{}\t{}\n", dst[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)], src[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)], face_gray_stencil[FaceCoordsCellGray::CELL_GRAY_C]);
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t i = 0; i < rowsize - 2; ++i)
  {
    for (size_t j = 0; j  < inner_rowsize - 2; ++j)
    {
      tmp = face_blue_stencil[FaceCoordsCellBlue::CELL_BLUE_C] * src[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)];

      out << fmt::format("{}\t{}\t{}\n", dst[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)], src[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)], face_blue_stencil[FaceCoordsCellBlue::CELL_BLUE_C]);
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, saveOperator_tmpl, saveOperator)

template< typename ValueType, size_t Level >
inline void printFunctionMemory(Face& face, const PrimitiveDataID<FaceBubbleFunctionMemory< ValueType >, Face> &dstId){
  using namespace std;
  ValueType* faceMemory = face.getData(dstId)->getPointer( Level );
  uint_t verticesPerDge = hhg::levelinfo::num_microvertices_per_edge(Level);
  cout << setfill('=') << setw(100) << "" << endl;
  cout << face << std::left << setprecision(1) << fixed << setfill(' ') << endl;
  std::cout << "Cell Blue: " << std::endl;
  for (size_t i = 0; i < verticesPerDge-2; ++i) {
    for (size_t j = 0; j < verticesPerDge-2 - i; ++j) {
      cout << setw(5) << faceMemory[FaceCoordsCellBlue::index<Level>(i, j, FaceCoordsCellBlue::CELL_BLUE_C)] << "|";
    }
    std::cout << std::endl;
  }
  cout << "Cell Gray: " << std::endl;
  for (size_t i = 0; i < verticesPerDge-1; ++i) {
    for (size_t j = 0; j < verticesPerDge-1 - i; ++j) {
      cout << setw(5) << faceMemory[FaceCoordsCellGray::index<Level>(i, j, FaceCoordsCellGray::CELL_GRAY_C)] << "|";
    }
    std::cout << std::endl;
  }
  cout << setw(100) << setfill(' ') << endl;

}

}// namespace P1BubbleFace
}// namespace hhg
