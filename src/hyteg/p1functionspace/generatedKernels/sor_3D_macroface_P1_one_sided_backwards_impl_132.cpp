
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "sor_3D_macroface_P1_one_sided_backwards_impl_132.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void sor_3D_macroface_P1_one_sided_backwards_impl_132_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int32_t level, double relax, std::map< hyteg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   for (int ctr_2 = (1 << (level)) - 1; ctr_2 >= 1; ctr_2 += -1)
   {
      // inner triangle
      for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 >= 1; ctr_1 += -1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_backwards_impl_132(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int32_t level, double relax, std::map< hyteg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_backwards_impl_132_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg