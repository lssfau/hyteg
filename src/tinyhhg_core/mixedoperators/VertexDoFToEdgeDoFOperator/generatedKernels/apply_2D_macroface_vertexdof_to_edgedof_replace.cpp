
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_2D_macroface_vertexdof_to_edgedof_replace.hpp"

namespace hhg {
namespace VertexDoFToEdgeDoF {
namespace generated {

static void apply_2D_macroface_vertexdof_to_edgedof_replace_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_vertexFaceSrc, double const * RESTRICT const _data_vertexToDiagonalEdgeFaceStencil, double const * RESTRICT const _data_vertexToHorizontalEdgeFaceStencil, double const * RESTRICT const _data_vertexToVerticalEdgeFaceStencil, int32_t level)
{
   const double xi_28 = _data_vertexToDiagonalEdgeFaceStencil[3];
   const double xi_29 = _data_vertexToDiagonalEdgeFaceStencil[2];
   const double xi_30 = _data_vertexToDiagonalEdgeFaceStencil[0];
   const double xi_31 = _data_vertexToDiagonalEdgeFaceStencil[1];
   const double xi_72 = _data_vertexToVerticalEdgeFaceStencil[3];
   const double xi_73 = _data_vertexToVerticalEdgeFaceStencil[0];
   const double xi_74 = _data_vertexToVerticalEdgeFaceStencil[2];
   const double xi_75 = _data_vertexToVerticalEdgeFaceStencil[1];
   const double xi_88 = _data_vertexToHorizontalEdgeFaceStencil[1];
   const double xi_89 = _data_vertexToHorizontalEdgeFaceStencil[3];
   const double xi_90 = _data_vertexToHorizontalEdgeFaceStencil[2];
   const double xi_91 = _data_vertexToHorizontalEdgeFaceStencil[0];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_35 = xi_28*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_36 = xi_29*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_37 = xi_30*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_38 = xi_31*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_35 + xi_36 + xi_37 + xi_38;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_79 = xi_28*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_80 = xi_29*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_81 = xi_30*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_82 = xi_31*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_83 = xi_72*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_84 = xi_73*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_85 = xi_74*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_86 = xi_75*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_79 + xi_80 + xi_81 + xi_82;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_83 + xi_84 + xi_85 + xi_86;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_47 = xi_72*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_48 = xi_73*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_49 = xi_74*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_50 = xi_75*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_47 + xi_48 + xi_49 + xi_50;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_99 = xi_88*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_89*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_101 = xi_90*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_102 = xi_91*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_103 = xi_28*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_104 = xi_29*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_105 = xi_30*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_106 = xi_31*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_99;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_103 + xi_104 + xi_105 + xi_106;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_15 = xi_88*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_16 = xi_89*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_17 = xi_90*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_18 = xi_91*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_19 = xi_28*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = xi_29*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_21 = xi_30*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_22 = xi_31*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_25 = xi_72*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_26 = xi_73*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_23 = xi_74*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_24 = xi_75*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_15 + xi_16 + xi_17 + xi_18;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_19 + xi_20 + xi_21 + xi_22;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25 + xi_26;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_119 = xi_88*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_120 = xi_89*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_121 = xi_90*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_122 = xi_91*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_72*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_124 = xi_73*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_125 = xi_74*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_126 = xi_75*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_119 + xi_120 + xi_121 + xi_122;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_123 + xi_124 + xi_125 + xi_126;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_59 = xi_88*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_89*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_61 = xi_90*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_62 = xi_91*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_59 + xi_60 + xi_61 + xi_62;
         }
      }
   }
}


void apply_2D_macroface_vertexdof_to_edgedof_replace(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_vertexFaceSrc, double const * RESTRICT const _data_vertexToDiagonalEdgeFaceStencil, double const * RESTRICT const _data_vertexToHorizontalEdgeFaceStencil, double const * RESTRICT const _data_vertexToVerticalEdgeFaceStencil, int32_t level)
{
    switch( level )
    {

    default:
        apply_2D_macroface_vertexdof_to_edgedof_replace_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_vertexFaceSrc, _data_vertexToDiagonalEdgeFaceStencil, _data_vertexToHorizontalEdgeFaceStencil, _data_vertexToVerticalEdgeFaceStencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace VertexDoFToEdgeDoF
} // namespace hhg