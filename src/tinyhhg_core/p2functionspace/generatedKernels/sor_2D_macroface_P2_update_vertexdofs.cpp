
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsP2MacroFace2D.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

static void sor_2D_macroface_P2_update_vertexdofs_level_2(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 6] - xi_10*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] - xi_13*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] - xi_15*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 6] - xi_16*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_17*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_18*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] - xi_6*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_7*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_8*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_9*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_3(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 10] - xi_10*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] - xi_13*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] - xi_15*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 10] - xi_16*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_17*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_18*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] - xi_6*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_7*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_8*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_9*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_4(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 18] - xi_10*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] - xi_13*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] - xi_15*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 18] - xi_16*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_17*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_18*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] - xi_6*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_7*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_8*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_9*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_5(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 34] - xi_10*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] - xi_13*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] - xi_15*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 34] - xi_16*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_17*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_18*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] - xi_6*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_7*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_8*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_9*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_6(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 66] - xi_10*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] - xi_13*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] - xi_15*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 66] - xi_16*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_17*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_18*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] - xi_6*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_7*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_8*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_9*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_7(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 130] - xi_10*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] - xi_13*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] - xi_15*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 130] - xi_16*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_17*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_18*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] - xi_6*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_7*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_8*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_9*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_8(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 258] - xi_10*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] - xi_13*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] - xi_15*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 258] - xi_16*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_17*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_18*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] - xi_6*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_7*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_8*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_9*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_9(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 514] - xi_10*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] - xi_13*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] - xi_15*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 514] - xi_16*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_17*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_18*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] - xi_6*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_7*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_8*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_9*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_10(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1026] - xi_10*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] - xi_13*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] - xi_15*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1026] - xi_16*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_17*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_18*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] - xi_6*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_7*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_8*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_9*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_11(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2050] - xi_10*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] - xi_13*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] - xi_15*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2050] - xi_16*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_17*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_18*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] - xi_6*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_7*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_8*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_9*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_12(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4098] - xi_10*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] - xi_13*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] - xi_15*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4098] - xi_16*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_17*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_18*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] - xi_6*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_7*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_8*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_9*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_13(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8194] - xi_10*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] - xi_13*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] - xi_15*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8194] - xi_16*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_17*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_18*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] - xi_6*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_7*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_8*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_9*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_14(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16386] - xi_10*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] - xi_13*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] - xi_15*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16386] - xi_16*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_17*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_18*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] - xi_6*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_7*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_8*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_9*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}

static void sor_2D_macroface_P2_update_vertexdofs_level_any(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, int64_t level, double relax)
{
   const double xi_0 = _data_vertex_stencil_at_vertex[3];
   const double xi_20 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_vertex[1];
   const double xi_2 = _data_edge_stencil_at_vertex[0];
   const double xi_3 = _data_edge_stencil_at_vertex[10];
   const double xi_4 = _data_edge_stencil_at_vertex[11];
   const double xi_5 = _data_edge_stencil_at_vertex[9];
   const double xi_6 = _data_edge_stencil_at_vertex[3];
   const double xi_7 = _data_edge_stencil_at_vertex[4];
   const double xi_8 = _data_edge_stencil_at_vertex[2];
   const double xi_9 = _data_edge_stencil_at_vertex[6];
   const double xi_10 = _data_edge_stencil_at_vertex[7];
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
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_20*(-xi_1*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1] - xi_10*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1] - xi_13*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_15*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_16*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_17*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_18*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1] - xi_4*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1] - xi_5*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_6*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_8*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_9*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_2D_macroface_P2_update_vertexdofs(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, int64_t level, double relax)
{
    switch( level )
    {
    case 2:
        sor_2D_macroface_P2_update_vertexdofs_level_2(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 3:
        sor_2D_macroface_P2_update_vertexdofs_level_3(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 4:
        sor_2D_macroface_P2_update_vertexdofs_level_4(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 5:
        sor_2D_macroface_P2_update_vertexdofs_level_5(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 6:
        sor_2D_macroface_P2_update_vertexdofs_level_6(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 7:
        sor_2D_macroface_P2_update_vertexdofs_level_7(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 8:
        sor_2D_macroface_P2_update_vertexdofs_level_8(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 9:
        sor_2D_macroface_P2_update_vertexdofs_level_9(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 10:
        sor_2D_macroface_P2_update_vertexdofs_level_10(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 11:
        sor_2D_macroface_P2_update_vertexdofs_level_11(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 12:
        sor_2D_macroface_P2_update_vertexdofs_level_12(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 13:
        sor_2D_macroface_P2_update_vertexdofs_level_13(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    case 14:
        sor_2D_macroface_P2_update_vertexdofs_level_14(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, relax);
        break;
    default:
        sor_2D_macroface_P2_update_vertexdofs_level_any(_data_edgeFaceDst, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, level, relax);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg