#pragma once

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/form/Form.hpp"
#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"

class p1_diffusion_cell_integral_0_otherwise;

namespace hhg {

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P1FenicsForm : public Form
{
 public:
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const override
   {
      Matrix3r localStiffnessMatrix;
      real_t   fenicsCoords[6];
      fenicsCoords[0] = coords[0][0];
      fenicsCoords[1] = coords[0][1];
      fenicsCoords[2] = coords[1][0];
      fenicsCoords[3] = coords[1][1];
      fenicsCoords[4] = coords[2][0];
      fenicsCoords[5] = coords[2][1];
      UFCOperator2D gen;
      gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
   }

   bool assemble2D() const { return !std::is_same< UFCOperator2D, hhg::fenics::NoAssemble >::value; }

   bool assemble3D() const { return !std::is_same< UFCOperator3D, hhg::fenics::NoAssemble >::value; }

   bool assembly2DDefined() const { return !std::is_same< UFCOperator2D, hhg::fenics::UndefinedAssembly >::value; }

   bool assembly3DDefined() const { return !std::is_same< UFCOperator3D, hhg::fenics::UndefinedAssembly >::value; }
};

} // namespace hhg