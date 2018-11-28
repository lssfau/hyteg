#pragma once

#include <functional>

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/types/matrix.hpp"

namespace hhg {
namespace fenics {

using walberla::real_c;

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

/// Use this UFCOperator to assemble the zero-matrix.
class NoAssemble {
 public:
  void tabulate_tensor(real_t *,
                       const real_t * const *,
                       const real_t *,
                       int) const { }
};

/// Use this UFCOperator to indicate that no assembly is defined at all.
class UndefinedAssembly {
public:
    void tabulate_tensor(real_t *,
                         const real_t * const *,
                         const real_t *,
                         int) const { WALBERLA_ABORT( "Assembly undefined." ); }
};

/// Dummy UFCOperator for a 10x10 (i.e. tet, P2) stiffness matrix
class Dummy10x10Assembly
{
public:

    Dummy10x10Assembly() : stiffnessMatrix_( real_c(0) ) {}
    Dummy10x10Assembly( const real_t & constant ) : stiffnessMatrix_( constant ) {}

    void tabulate_tensor(real_t * A,
                         const real_t * const *,
                         const real_t *,
                         int) const
    {
      for ( uint_t i = 0; i < 100; i++ )
        A[i] = stiffnessMatrix_.data()[i];
    }

private:

    Matrix10r stiffnessMatrix_;
};


typedef std::function<void(real_t *,
                      const real_t * const *,
                      const real_t *,
                      int cell_orientation)> TabulateTensor;


}
}
