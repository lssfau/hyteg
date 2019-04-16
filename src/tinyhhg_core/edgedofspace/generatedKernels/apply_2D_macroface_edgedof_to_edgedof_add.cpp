
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsEdgeToEdgeMacroFace2D.hpp"

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

static void apply_2D_macroface_edgedof_to_edgedof_add_level_2(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         _data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_3(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         _data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_4(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         _data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_5(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         _data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_6(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         _data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_7(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         _data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_8(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         _data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_9(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         _data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_10(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         _data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_11(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         _data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_12(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         _data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_13(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         _data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_14(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         _data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}

static void apply_2D_macroface_edgedof_to_edgedof_add_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil, int64_t level)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_138 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_139 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_140 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_141 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_142 = _data_edgeToHorizontalEdgeFaceStencil[2];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
      {
         const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
         _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
      }
      // bottom right vertex
      for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
      {
         const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
   for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_156 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_162 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_151 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_152 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_153 = xi_140*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_154 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_155 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_157 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_158 = xi_38*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_159 = xi_39*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_160 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_161 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
         _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
      {
         const double xi_23 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_140*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_21 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_22 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
         _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
         _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_130 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_125 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_126 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_127 = xi_140*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_128 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_129 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_131 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_133 = xi_93*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_134 = xi_94*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
         _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
      }
   }
   for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_80 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_75 = xi_138*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_76 = xi_139*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_77 = xi_140*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_78 = xi_141*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_79 = xi_142*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
      }
   }
   {
      
   }
}


void apply_2D_macroface_edgedof_to_edgedof_add(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil, int64_t level)
{
    switch( level )
    {
    case 2:
        apply_2D_macroface_edgedof_to_edgedof_add_level_2(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 3:
        apply_2D_macroface_edgedof_to_edgedof_add_level_3(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 4:
        apply_2D_macroface_edgedof_to_edgedof_add_level_4(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 5:
        apply_2D_macroface_edgedof_to_edgedof_add_level_5(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 6:
        apply_2D_macroface_edgedof_to_edgedof_add_level_6(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 7:
        apply_2D_macroface_edgedof_to_edgedof_add_level_7(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 8:
        apply_2D_macroface_edgedof_to_edgedof_add_level_8(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 9:
        apply_2D_macroface_edgedof_to_edgedof_add_level_9(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 10:
        apply_2D_macroface_edgedof_to_edgedof_add_level_10(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 11:
        apply_2D_macroface_edgedof_to_edgedof_add_level_11(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 12:
        apply_2D_macroface_edgedof_to_edgedof_add_level_12(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 13:
        apply_2D_macroface_edgedof_to_edgedof_add_level_13(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    case 14:
        apply_2D_macroface_edgedof_to_edgedof_add_level_14(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil);
        break;
    default:
        apply_2D_macroface_edgedof_to_edgedof_add_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg