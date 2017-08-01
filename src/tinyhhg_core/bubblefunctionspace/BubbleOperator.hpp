#pragma once

#include <fmt/format.h>
#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "BubbleDataHandling.hpp"

#include "tinyhhg_core/fenics.hpp"

#include "generated/bubble_diffusion.h"

#include "BubbleMemory.hpp"
#include "BubbleFace.hpp"

namespace hhg
{

template<class UFCOperator>
class BubbleOperator : public Operator< BubbleFunction, BubbleFunction >
{
public:
  BubbleOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel)
  {
    FaceBubbleStencilMemoryDataHandling faceBubbleStencilMemoryDataHandling(minLevel_, maxLevel_);
    faceStencilID_ = storage->addFaceData(faceBubbleStencilMemoryDataHandling, "BubbleOperatorFaceStencil");

    for (uint_t level = minLevel_; level <= maxLevel_; ++level)
    {

      for (auto& it : storage_->getFaces()) {
        Face& face = *it.second;

        auto& face_stencil_stack = face.getData(faceStencilID_)->data[level];

        real_t local_stiffness_gray[1][1];
        real_t local_stiffness_blue[1][1];
        compute_local_stiffness(face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(face, level, local_stiffness_blue, fenics::BLUE);

        face_stencil_stack[fenics::GRAY][0] = local_stiffness_gray[0][0];
        face_stencil_stack[fenics::BLUE][0] = local_stiffness_blue[0][0];
      }
    }
  }

  ~BubbleOperator()
  {
  }

private:

  void apply_impl(BubbleFunction& src, BubbleFunction& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        BubbleFace::apply(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }
  }

  void smooth_gs_impl(BubbleFunction& dst, BubbleFunction& rhs, size_t level, DoFType flag)
  {
    WALBERLA_ASSERT(false, "BubbleOperator::smooth_gs is not implemented!");
  }

private:
  PrimitiveDataID<FaceBubbleStencilMemory, Face> faceStencilID_;

  void compute_local_stiffness(const Face &face, size_t level, real_t local_stiffness[1][1], fenics::ElementType element_type) {
    real_t A[1];
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator gen;
    gen.tabulate_tensor(A, NULL, coords, 0);

    local_stiffness[0][0] = A[0];
  }

};

typedef BubbleOperator<bubble_diffusion_cell_integral_0_otherwise> BubbleLaplaceOperator;
}
