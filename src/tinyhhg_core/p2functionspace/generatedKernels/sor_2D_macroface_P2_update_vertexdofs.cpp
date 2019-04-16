
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsP2MacroFace2D.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

static void sor_2D_macroface_P2_update_vertexdofs_level_2(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 4; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 6];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 6];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_3(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 10];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 10];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_4(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 16; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 18];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 18];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_5(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 32; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 34];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 34];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_6(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 64; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 66];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 66];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_7(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 128; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 130];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 130];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_8(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 256; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 258];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 258];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_9(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 512; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 514];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 514];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_10(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 1024; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1026];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1026];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_11(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 2048; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2050];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2050];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_12(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 4096; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4098];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4098];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_13(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 8192; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8194];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8194];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_14(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < 16384; ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16386];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16386];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, int64_t level, double relax)
{
   const double xi_23 = 1.0;
   const double xi_24 = -relax;
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_21 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[10];
   const double xi_3 = _data_edge_stencil_at_vertex[4];
   const double xi_4 = _data_edge_stencil_at_vertex[7];
   const double xi_5 = _data_edge_stencil_at_vertex[0];
   const double xi_6 = _data_edge_stencil_at_vertex[9];
   const double xi_7 = _data_edge_stencil_at_vertex[3];
   const double xi_8 = _data_edge_stencil_at_vertex[6];
   const double xi_9 = _data_edge_stencil_at_vertex[11];
   const double xi_10 = _data_edge_stencil_at_vertex[2];
   const double xi_11 = _data_edge_stencil_at_vertex[8];
   const double xi_12 = _data_edge_stencil_at_vertex[5];
   const double xi_13 = _data_vertex_stencil_at_vertex[2];
   const double xi_14 = _data_vertex_stencil_at_vertex[5];
   const double xi_15 = _data_vertex_stencil_at_vertex[0];
   const double xi_16 = _data_vertex_stencil_at_vertex[6];
   const double xi_17 = _data_vertex_stencil_at_vertex[1];
   const double xi_18 = _data_vertex_stencil_at_vertex[4];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_43 = _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = -xi_1*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_35 = -xi_2*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = -xi_3*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_37 = -xi_4*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_5*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_39 = -xi_6*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_40 = -xi_7*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_41 = -xi_8*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = -xi_9*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = -xi_10*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_27 = -xi_11*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_28 = -xi_12*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_29 = -xi_13*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_30 = -xi_14*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_31 = -xi_15*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_32 = -xi_16*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_33 = -xi_17*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_34 = -xi_18*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_21*(xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_23 + xi_24)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_2D_macroface_P2_update_vertexdofs(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * const _data_vertex_stencil_at_vertex, int64_t level, double relax)
{
    switch( level )
    {
    case 2:
        sor_2D_macroface_P2_update_vertexdofs_level_2(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 3:
        sor_2D_macroface_P2_update_vertexdofs_level_3(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 4:
        sor_2D_macroface_P2_update_vertexdofs_level_4(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 5:
        sor_2D_macroface_P2_update_vertexdofs_level_5(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 6:
        sor_2D_macroface_P2_update_vertexdofs_level_6(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 7:
        sor_2D_macroface_P2_update_vertexdofs_level_7(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 8:
        sor_2D_macroface_P2_update_vertexdofs_level_8(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 9:
        sor_2D_macroface_P2_update_vertexdofs_level_9(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 10:
        sor_2D_macroface_P2_update_vertexdofs_level_10(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 11:
        sor_2D_macroface_P2_update_vertexdofs_level_11(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 12:
        sor_2D_macroface_P2_update_vertexdofs_level_12(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 13:
        sor_2D_macroface_P2_update_vertexdofs_level_13(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 14:
        sor_2D_macroface_P2_update_vertexdofs_level_14(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    default:
        sor_2D_macroface_P2_update_vertexdofs_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, level, relax);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg