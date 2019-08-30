
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_2D_macroface_edgedof_to_vertexdof_replace.hpp"

namespace hyteg {
namespace EdgeDoFToVertexDoF {
namespace generated {

static void apply_2D_macroface_edgedof_to_vertexdof_replace_level_any(double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeToVertexFaceStencil, double * RESTRICT _data_p1FaceDst, int32_t level)
{
   const double xi_0 = _data_edgeToVertexFaceStencil[1];
   const double xi_1 = _data_edgeToVertexFaceStencil[10];
   const double xi_2 = _data_edgeToVertexFaceStencil[4];
   const double xi_3 = _data_edgeToVertexFaceStencil[7];
   const double xi_4 = _data_edgeToVertexFaceStencil[0];
   const double xi_5 = _data_edgeToVertexFaceStencil[9];
   const double xi_6 = _data_edgeToVertexFaceStencil[3];
   const double xi_7 = _data_edgeToVertexFaceStencil[6];
   const double xi_8 = _data_edgeToVertexFaceStencil[11];
   const double xi_9 = _data_edgeToVertexFaceStencil[2];
   const double xi_10 = _data_edgeToVertexFaceStencil[8];
   const double xi_11 = _data_edgeToVertexFaceStencil[5];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_15 = xi_0*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_16 = xi_1*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_19 = xi_2*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_20 = xi_3*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_4*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = xi_5*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_6*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_24 = xi_7*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_8*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = xi_9*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_17 = xi_10*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_11*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26;
      }
   }
}


void apply_2D_macroface_edgedof_to_vertexdof_replace(double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeToVertexFaceStencil, double * RESTRICT _data_p1FaceDst, int32_t level)
{
    switch( level )
    {

    default:
        apply_2D_macroface_edgedof_to_vertexdof_replace_level_any(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToVertexFaceStencil, _data_p1FaceDst, level);
        break;
    }
}
    

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hyteg