#pragma once

#include <fmt/format.h>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/operator.hpp"

#include "BubbleDataHandling.hpp"

#include "generated/bubble_diffusion.h"

#include "BubbleMemory.hpp"
#include "BubbleFace.hpp"

namespace hhg
{

namespace BubbleSpace {
enum ElementType {
  GRAY,
  BLUE
};

void compute_micro_coords(const Face &face, size_t level, real_t coords[6], ElementType element_type) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  Point3D d0 = (face.coords[1] - face.coords[0]) / walberla::real_c((rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0]) / walberla::real_c((rowsize - 1));

  real_t orientation = 1.0;

  if (element_type == BLUE) {
    orientation = -1.0;
  }

  coords[0] = 0.0;
  coords[1] = 0.0;
  coords[2] = orientation * d0[0];
  coords[3] = orientation * d0[1];
  coords[4] = orientation * d2[0];
  coords[5] = orientation * d2[1];
}

template<class UFCOperator>
void compute_local_stiffness(const Face &face, size_t level, real_t local_stiffness[1][1], ElementType element_type) {
  real_t A[1];
  real_t coords[6];
  compute_micro_coords(face, level, coords, element_type);
  UFCOperator gen;
  gen.tabulate_tensor(A, NULL, coords, 0);

  local_stiffness[0][0] = A[0];
}
}

template<class UFCOperator>
class BubbleOperator : public Operator
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
        BubbleSpace::compute_local_stiffness<UFCOperator>(face, level, local_stiffness_gray, BubbleSpace::GRAY);
        BubbleSpace::compute_local_stiffness<UFCOperator>(face, level, local_stiffness_blue, BubbleSpace::BLUE);

        face_stencil_stack[BubbleSpace::GRAY][0] = local_stiffness_gray[0][0];
        face_stencil_stack[BubbleSpace::BLUE][0] = local_stiffness_blue[0][0];
      }
    }
  }

  ~BubbleOperator()
  {
  }

  void apply(BubbleFunction& src, BubbleFunction& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        BubbleFace::apply(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }
  }

  void smooth_gs(P1Function& dst, P1Function& rhs, size_t level, DoFType flag)
  {
    WALBERLA_ASSERT(false, "BubbleOperator::smooth_gs is not implemented!");
  }

 private:
  PrimitiveDataID<FaceBubbleStencilMemory, Face> faceStencilID_;

};

typedef BubbleOperator<bubble_diffusion_cell_integral_0_otherwise> BubbleLaplaceOperator;
}
