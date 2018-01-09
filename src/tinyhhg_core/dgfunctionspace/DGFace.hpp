
#pragma once

#include "DGFaceIndex.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"

namespace hhg {
namespace DGFace {

template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Face &face,
                          const PrimitiveDataID < FunctionMemory< ValueType >, Face> &dstId,
                          uint_t& num)
{
  ValueType *dst = face.getData(dstId)->getPointer(Level);
  //the outermost gray cells belong to the edge
  const uint_t grayRowsize = levelinfo::num_microvertices_per_edge(Level) - 1;
  uint_t inner_rowsize = grayRowsize - 2;
  for(uint_t col = 1; col < grayRowsize - 2; ++col){
    for(uint_t row = 1; row < inner_rowsize; ++row){
      dst[indexDGFaceFromVertex< Level >(col,row,stencilDirection::CELL_GRAY_NE)] = num;
      ++num;
    }
    --inner_rowsize;
  }

  const uint_t blueRowsize = levelinfo::num_microvertices_per_edge(Level) - 2;
  inner_rowsize = blueRowsize;
  for(uint_t col = 0; col < blueRowsize; ++col){
    for(uint_t row = 0; row < inner_rowsize; ++row){
      dst[indexDGFaceFromGrayDGface< Level >(col,row,stencilDirection::CELL_BLUE_E)] = num;
      ++num;
    }
    --inner_rowsize;
  }

}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate )

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Face &face,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face>& faceMemoryId,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>>& srcMemoryIds,
                            std::function<ValueType(const hhg::Point3D &, const std::vector<ValueType>& f)> &expr) {

  auto faceMemory = face.getData(faceMemoryId)->getPointer( Level );

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0;

  x0 = face.coords[0];

  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));

  uint_t inner_rowsize = rowsize;

  std::vector<ValueType*> srcPtr;
  for(auto src : srcMemoryIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }

  std::vector<ValueType> srcVector(srcMemoryIds.size());

  // gray cells
  for (size_t j = 1; j < rowsize - 2; ++j) {

    x = x0 + 1.0/3.0 * (d0 + d2) + walberla::real_c(j) * d2 + d0;

    for (size_t i = 1; i < inner_rowsize - 3; ++i) {

      for (size_t k = 0; k < srcPtr.size(); ++k) {
        srcVector[k] = srcPtr[k][BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)];
      }

      faceMemory[BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] = expr(x, srcVector);
      x += d0;
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  // blue cells
  for (size_t j = 0; j < rowsize - 2; ++j) {

    x = x0 + 2.0/3.0 * (d0 + d2) + walberla::real_c(j) * d2;

    for (size_t i = 0; i < inner_rowsize - 2; ++i) {
      for (size_t k = 0; k < srcPtr.size(); ++k) {
        srcVector[k] = srcPtr[k][BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)];
      }

      faceMemory[BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] = expr(x, srcVector);
      x += d0;
    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, interpolateTmpl, interpolate)

template< typename ValueType, uint_t Level >
inline void addTmpl(Face &face,
                            const std::vector<ValueType>& scalars,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &srcIds,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto dst = face.getData(dstId)->getPointer(Level);

  // gray cells
  for (size_t j = 1; j < rowsize - 2; ++j) {
    for (size_t i = 1; i < inner_rowsize - 3; ++i) {

      auto cellIndex = BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C);

      ValueType tmp = 0.0;

      for (uint_t k = 0; k < srcIds.size(); ++k) {
        tmp += scalars[k]*face.getData(srcIds[k])->getPointer( Level )[cellIndex];
      }

      dst[cellIndex] += tmp;

    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  // blue cells
  for (size_t j = 0; j < rowsize - 2; ++j) {
    for (size_t i = 0; i < inner_rowsize - 2; ++i) {

      auto cellIndex = BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C);

      ValueType tmp = 0.0;

      for (uint_t k = 0; k < srcIds.size(); ++k) {
        tmp += scalars[k]*face.getData(srcIds[k])->getPointer( Level )[cellIndex];
      }

      dst[cellIndex] += tmp;

    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, addTmpl, add)

template< typename ValueType, uint_t Level >
inline void assignTmpl(Face &face,
                    const std::vector<ValueType>& scalars,
                    const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &srcIds,
                    const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId) {

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  auto dst = face.getData(dstId)->getPointer(Level);
  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }
  // gray cells
  for (size_t j = 1; j < rowsize - 2; ++j) {
    for (size_t i = 1; i < inner_rowsize - 3; ++i) {

      auto cellIndex = BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C);

      ValueType tmp = scalars[0]*srcPtr[0][cellIndex];

      for (uint_t k = 1; k < srcIds.size(); ++k) {
        tmp += scalars[k]*srcPtr[k][cellIndex];
      }

      dst[cellIndex] = tmp;

    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  // blue cells
  for (size_t j = 0; j < rowsize - 2; ++j) {
    for (size_t i = 0; i < inner_rowsize - 2; ++i) {

      auto cellIndex = BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C);

      ValueType tmp = scalars[0]*srcPtr[0][cellIndex];

      for (uint_t k = 1; k < srcIds.size(); ++k) {
        tmp += scalars[k]*srcPtr[k][cellIndex];
      }

      dst[cellIndex] = tmp;

    }
    --inner_rowsize;
  }
}

SPECIALIZE_WITH_VALUETYPE(void, assignTmpl, assign)

template< typename ValueType, uint_t Level >
inline void upwindTmpl(Face &face,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       const PrimitiveDataID < FunctionMemory< ValueType >, Face> &srcId,
                       const PrimitiveDataID < FunctionMemory< ValueType >, Face> &dstId,
                       const std::array<PrimitiveDataID< FunctionMemory< ValueType >, Face>, 2> &velocityIds,
                       UpdateType updateType) {

  using namespace vertexdof::macroface;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  // get memories
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );
  auto u = face.getData(velocityIds[0])->getPointer(Level);
  auto v = face.getData(velocityIds[1])->getPointer(Level);

  // get edge directions
  auto d0 = (face.coords[1] - face.coords[0]) / walberla::real_c(rowsize-1);
  auto d1 = (face.coords[2] - face.coords[1]) / walberla::real_c(rowsize-1);
  auto d2 = (face.coords[0] - face.coords[2]) / walberla::real_c(rowsize-1);

  // compute edge lengths
  real_t d0Length = d0.norm();
  real_t d1Length = d1.norm();
  real_t d2Length = d2.norm();

  // compute normals
  auto n_0 = d0.normal2D() / d0Length;
  auto n_1 = d1.normal2D() / d1Length;
  auto n_2 = d2.normal2D() / d2Length;

  real_t faceOrientation = math::faceOrientation2D(face.coords[0], face.coords[1], face.coords[2]);

  // correct normals if all normals point in wrong direction
  n_0 *= faceOrientation;
  n_1 *= faceOrientation;
  n_2 *= faceOrientation;

  real_t faceArea = std::pow(4.0, -walberla::real_c(Level)) * face.area;
  real_t faceAreaInv = 1.0 / faceArea;

  ValueType tmp;

  Point2D u_0, u_1, u_2;
  real_t un_0, un_1, un_2;
  real_t c_up_0, c_up_1, c_up_2;

  for (size_t j = 1; j < rowsize - 2; ++j)
  {
    for (size_t i = 1; i  < inner_rowsize - 3; ++i)
    {
      // evalate velocities
      u_0[0] = 0.5 * (u[indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] + u[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)]);
      u_0[1] = 0.5 * (v[indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] + v[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)]);

      u_1[0] = 0.5 * (u[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)] + u[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)]);
      u_1[1] = 0.5 * (v[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)] + v[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)]);

      u_2[0] = 0.5 * (u[indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] + u[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)]);
      u_2[1] = 0.5 * (v[indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] + v[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)]);

      un_0 = d0Length * u_0.dot(n_0);
      un_1 = d1Length * u_1.dot(n_1);
      un_2 = d2Length * u_2.dot(n_2);

      if (un_0 >= 0) {
        c_up_0 = src[BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)];
      }
      else {
        c_up_0 = src[indexDGFaceFromGrayDGface<Level>(i, j, stencilDirection::CELL_BLUE_S)];
      }

      if (un_1 >= 0) {
        c_up_1 = src[BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)];
      }
      else {
        c_up_1 = src[indexDGFaceFromGrayDGface<Level>(i, j, stencilDirection::CELL_BLUE_E)];
      }

      if (un_2 >= 0) {
        c_up_2 = src[BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)];
      }
      else {
        c_up_2 = src[indexDGFaceFromGrayDGface<Level>(i, j, stencilDirection::CELL_BLUE_W)];
      }

      tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
      tmp *= faceAreaInv;

      if (updateType == Replace) {
        dst[BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] = tmp;
      } else if (updateType == Add) {
        dst[BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] += tmp;
      }
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  // flip normals
  n_0 *= -1.0;
  n_1 *= -1.0;
  n_2 *= -1.0;

  for (size_t j = 0; j < rowsize - 2; ++j)
  {
    for (size_t i = 0; i  < inner_rowsize - 2; ++i)
    {
      // evalate velocities
      u_0[0] = 0.5 * (u[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)] + u[indexFromVertex<Level>(i+1, j+1, stencilDirection::VERTEX_C)]);
      u_0[1] = 0.5 * (v[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)] + v[indexFromVertex<Level>(i+1, j+1, stencilDirection::VERTEX_C)]);

      u_1[0] = 0.5 * (u[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)] + u[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)]);
      u_1[1] = 0.5 * (v[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)] + v[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)]);

      u_2[0] = 0.5 * (u[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)] + u[indexFromVertex<Level>(i+1, j+1, stencilDirection::VERTEX_C)]);
      u_2[1] = 0.5 * (v[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)] + v[indexFromVertex<Level>(i+1, j+1, stencilDirection::VERTEX_C)]);

      un_0 = d0Length * u_0.dot(n_0);
      un_1 = d1Length * u_1.dot(n_1);
      un_2 = d2Length * u_2.dot(n_2);

      if (un_0 >= 0) {
        c_up_0 = src[BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)];
      }
      else {
        c_up_0 = src[indexDGFaceFromBlueDGface<Level>(i, j, stencilDirection::CELL_GRAY_N)];
      }

      if (un_1 >= 0) {
        c_up_1 = src[BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)];
      }
      else {
        c_up_1 = src[indexDGFaceFromBlueDGface<Level>(i, j, stencilDirection::CELL_GRAY_W)];
      }

      if (un_2 >= 0) {
        c_up_2 = src[BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)];
      }
      else {
        c_up_2 = src[indexDGFaceFromBlueDGface<Level>(i, j, stencilDirection::CELL_GRAY_E)];
      }

      tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
      tmp *= faceAreaInv;

      if (updateType == Replace) {
        dst[BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] = tmp;
      } else if (updateType == Add) {
        dst[BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] += tmp;
      }
    }
    --inner_rowsize;
  }

}

SPECIALIZE_WITH_VALUETYPE( void, upwindTmpl, upwind )

template< typename ValueType, uint_t Level >
inline void projectP1Tmpl(Face &face,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       const PrimitiveDataID < FunctionMemory< ValueType >, Face> &srcId,
                       const PrimitiveDataID < FunctionMemory< ValueType >, Face> &dstId,
                       UpdateType updateType) {

  using namespace vertexdof::macroface;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  // get memories
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );

  ValueType tmp;

  for (size_t j = 1; j < rowsize - 2; ++j)
  {
    for (size_t i = 1; i  < inner_rowsize - 3; ++i)
    {
      // evalate velocities
      tmp = 1.0/3.0 * (src[indexFromVertex<Level>(i, j, stencilDirection::VERTEX_C)] + src[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)] + src[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)]);

      if (updateType == Replace) {
        dst[BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] = tmp;
      } else if (updateType == Add) {
        dst[BubbleFace::indexFaceFromGrayFace<Level>(i, j, stencilDirection::CELL_GRAY_C)] += tmp;
      }
    }
    --inner_rowsize;
  }

  inner_rowsize = rowsize;

  for (size_t j = 0; j < rowsize - 2; ++j)
  {
    for (size_t i = 0; i  < inner_rowsize - 2; ++i)
    {
      // evalate velocities
      tmp = 1.0/3.0 * (src[indexFromVertex<Level>(i, j+1, stencilDirection::VERTEX_C)] + src[indexFromVertex<Level>(i+1, j+1, stencilDirection::VERTEX_C)] + src[indexFromVertex<Level>(i+1, j, stencilDirection::VERTEX_C)]);

      if (updateType == Replace) {
        dst[BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] = tmp;
      } else if (updateType == Add) {
        dst[BubbleFace::indexFaceFromBlueFace<Level>(i, j, stencilDirection::CELL_BLUE_C)] += tmp;
      }
    }
    --inner_rowsize;
  }

}

SPECIALIZE_WITH_VALUETYPE( void, projectP1Tmpl, projectP1 )

}//namespace DGFace
}//namespace hhg
