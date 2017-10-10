#pragma once

#include <fmt/format.h>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/Operator.hpp"
#include "P1ToBubbleMemory.hpp"
#include "P1ToBubbleDataHandling.hpp"

#include "tinyhhg_core/bubblefunctionspace/BubbleFaceIndex.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleEdgeIndex.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics.hpp"

#include "generated/p1_to_bubble_div.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "P1ToBubbleFace.hpp"

namespace hhg
{

template<class UFCOperator>
class P1ToBubbleOperator : public Operator<P1Function< real_t >, BubbleFunction< real_t > >
{
 public:
  P1ToBubbleOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel)
  {
    auto faceP1ToBubbleStencilMemoryDataHandling = std::make_shared< FaceP1ToBubbleStencilMemoryDataHandling >(minLevel_, maxLevel_);

    storage->addFaceData(faceStencilID_, faceP1ToBubbleStencilMemoryDataHandling, "P1ToBubbleOperatorFaceStencil");

    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {

      // assemble face stencil
      for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        auto& face_stencil_stack = face.getData(faceStencilID_)->data[level];

        auto& face_gray_stencil = face_stencil_stack[0];
        auto& face_blue_stencil = face_stencil_stack[1];

        real_t local_stiffness_gray[1][3];
        real_t local_stiffness_blue[1][3];
        compute_local_stiffness(face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(face, level, local_stiffness_blue, fenics::BLUE);

        face_gray_stencil[P1Face::FaceCoordsCellGray::VERTEX_SW] = local_stiffness_gray[0][0];
        face_gray_stencil[P1Face::FaceCoordsCellGray::VERTEX_SE] = local_stiffness_gray[0][1];
        face_gray_stencil[P1Face::FaceCoordsCellGray::VERTEX_NW] = local_stiffness_gray[0][2];

        face_blue_stencil[P1Face::FaceCoordsCellBlue::VERTEX_SE] = local_stiffness_blue[0][2];
        face_blue_stencil[P1Face::FaceCoordsCellBlue::VERTEX_NW] = local_stiffness_blue[0][1];
        face_blue_stencil[P1Face::FaceCoordsCellBlue::VERTEX_NE] = local_stiffness_blue[0][0];
      }
    }
  }

  ~P1ToBubbleOperator()
  {
  }

  void apply_impl(P1Function< real_t > & src, BubbleFunction< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        P1ToBubbleFace::apply(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }
  }

  const PrimitiveDataID<FaceP1ToBubbleStencilMemory, Face> &getFaceStencilID() const { return faceStencilID_; }

 private:
  PrimitiveDataID<FaceP1ToBubbleStencilMemory, Face> faceStencilID_;

  void compute_local_stiffness(const Face &face, size_t level, real_t local_stiffness[1][3], fenics::ElementType element_type) {
    real_t A[3];
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(A, NULL, coords, 0);

    for (size_t j = 0; j < 3; ++j) {
      size_t i = 0;
      local_stiffness[i][j] = A[j];
    }
  }
};

typedef P1ToBubbleOperator<p1_to_bubble_div_cell_integral_0_otherwise> P1ToBubbleDivxOperator;
typedef P1ToBubbleOperator<p1_to_bubble_div_cell_integral_1_otherwise> P1ToBubbleDivyOperator;

}
