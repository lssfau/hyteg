#pragma once

#include <functional>

namespace hhg
{
namespace fenics {

enum ElementType {
  GRAY,
  BLUE
};

inline void compute_micro_coords(const Face &face, size_t level, real_t coords[6], ElementType element_type) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(level);
  Point3D d0 = (face.coords[1] - face.coords[0])/walberla::real_c((rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/walberla::real_c((rowsize - 1));

  real_t orientation = 1.0;

  if (element_type==BLUE) {
    orientation = -1.0;
  }

  coords[0] = 0.0;
  coords[1] = 0.0;
  coords[2] = orientation*d0[0];
  coords[3] = orientation*d0[1];
  coords[4] = orientation*d2[0];
  coords[5] = orientation*d2[1];
}

class NoAssemble {
 public:
  void tabulate_tensor(real_t * A,
                       const real_t * const * w,
                       const real_t * coordinate_dofs,
                       int cell_orientation) const { }
};

typedef std::function<void(real_t *,
                      const real_t * const *,
                      const real_t *,
                      int cell_orientation)> TabulateTensor;


}
}
