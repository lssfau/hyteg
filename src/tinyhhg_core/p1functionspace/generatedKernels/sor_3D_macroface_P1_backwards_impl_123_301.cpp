
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "sor_3D_macroface_P1_backwards_impl_123_301.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void sor_3D_macroface_P1_backwards_impl_123_301_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double * RESTRICT _data_vertexFaceDst_gl1, double const * RESTRICT const _data_vertexFaceRhs, int32_t level, double relax, std::map< hyteg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0, std::map< hyteg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_1)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_2 = v2v_cell_stencil_fused_face_1[{ 0, 0, 0 }];
   const double xi_24 = 1 / (xi_1 + xi_2);
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_12 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_13 = v2v_cell_stencil_fused_face_1[{ -1, 0, 0 }];
   const double xi_14 = v2v_cell_stencil_fused_face_1[{ -1, 0, 1 }];
   const double xi_15 = v2v_cell_stencil_fused_face_1[{ -1, 1, -1 }];
   const double xi_16 = v2v_cell_stencil_fused_face_1[{ -1, 1, 0 }];
   const double xi_17 = v2v_cell_stencil_fused_face_1[{ 0, 0, -1 }];
   const double xi_18 = v2v_cell_stencil_fused_face_1[{ 0, 0, 1 }];
   const double xi_19 = v2v_cell_stencil_fused_face_1[{ 0, 1, -1 }];
   const double xi_20 = v2v_cell_stencil_fused_face_1[{ 0, 1, 0 }];
   const double xi_21 = v2v_cell_stencil_fused_face_1[{ 1, 0, -1 }];
   const double xi_22 = v2v_cell_stencil_fused_face_1[{ 1, 0, 0 }];
   for (int ctr_2 = (1 << (level)) - 1; ctr_2 >= 1; ctr_2 += -1)
   {
      // inner triangle
      for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 >= 1; ctr_1 += -1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_24*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_11*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_13*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_15*_data_vertexFaceDst_gl1[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_16*_data_vertexFaceDst_gl1[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_17*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_18*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_19*_data_vertexFaceDst_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_20*_data_vertexFaceDst_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_21*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_22*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_3*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_5*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_6*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_7*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_backwards_impl_123_301(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double * RESTRICT _data_vertexFaceDst_gl1, double const * RESTRICT const _data_vertexFaceRhs, int32_t level, double relax, std::map< hyteg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0, std::map< hyteg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_1)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_backwards_impl_123_301_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg