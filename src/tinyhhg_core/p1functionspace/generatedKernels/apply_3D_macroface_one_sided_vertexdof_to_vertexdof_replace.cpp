
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsVertexToVertexMacroFace3D.hpp"

namespace hhg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_012(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_3 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_4 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_5 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_6 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_7 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_8 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_9 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_10 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_013(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_3 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_4 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_5 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_6 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_7 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_8 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_9 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_10 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_021(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_3 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_4 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_5 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_6 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_7 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_8 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_9 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_10 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_023(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_3 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_4 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_5 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_6 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_7 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_8 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_9 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_10 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_031(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_3 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_4 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_5 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_6 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_7 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_8 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_9 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_10 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_032(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_3 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_4 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_5 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_6 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_7 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_8 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_9 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_10 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_102(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_3 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_4 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_5 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_6 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_7 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_8 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_9 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_10 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_103(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_3 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_4 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_5 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_6 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_7 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_8 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_9 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_10 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_120(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_3 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_4 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_5 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_6 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_7 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_8 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_9 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_10 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_123(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_3 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_4 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_5 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_6 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_7 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_8 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_9 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_10 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_130(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_3 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_4 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_5 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_6 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_7 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_8 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_9 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_10 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_132(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_3 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_4 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_5 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_6 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_7 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_8 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_9 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_10 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_201(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_3 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_4 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_5 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_6 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_7 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_8 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_9 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_10 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_203(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_3 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_4 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_5 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_6 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_7 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_8 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_9 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_10 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_210(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_3 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_4 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_5 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_6 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_7 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_8 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_9 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_10 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_213(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_3 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_4 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_5 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_6 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_7 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_8 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_9 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_10 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_230(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_3 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_4 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_5 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_6 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_7 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_8 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_9 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_10 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_231(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_3 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_4 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_5 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_6 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_7 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_8 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_9 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_10 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_301(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_3 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_4 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_5 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_6 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_7 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_8 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_9 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_10 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_302(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_3 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_4 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_5 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_6 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_7 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_8 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_9 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_10 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_310(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_3 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_4 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_5 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_6 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_7 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_8 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_9 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_10 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_312(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_3 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_4 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_5 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_6 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_7 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_8 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_9 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_10 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_320(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ 1, -1, 1 }];
   const double xi_2 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_3 = p1FaceStencil[{ 1, 0, 0 }];
   const double xi_4 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_5 = p1FaceStencil[{ 0, 1, 0 }];
   const double xi_6 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_7 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_8 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_9 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_10 = p1FaceStencil[{ 0, 0, 1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

static void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_321(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[{ -1, 1, -1 }];
   const double xi_2 = p1FaceStencil[{ 0, -1, 0 }];
   const double xi_3 = p1FaceStencil[{ -1, 0, 0 }];
   const double xi_4 = p1FaceStencil[{ 0, 0, -1 }];
   const double xi_5 = p1FaceStencil[{ -1, 1, 0 }];
   const double xi_6 = p1FaceStencil[{ 0, 1, -1 }];
   const double xi_7 = p1FaceStencil[{ 1, -1, 0 }];
   const double xi_8 = p1FaceStencil[{ 0, -1, 1 }];
   const double xi_9 = p1FaceStencil[{ 1, 0, -1 }];
   const double xi_10 = p1FaceStencil[{ -1, 0, 1 }];
   const double xi_11 = p1FaceStencil[{ 0, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_14 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_17 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_18 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_19 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_21 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_22 = xi_7*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_23 = xi_8*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = xi_9*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_15 = xi_10*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = xi_11*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24;
      }
   }
}

void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int64_t level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2, std::map< hhg::indexing::IndexIncrement, double > p1FaceStencil)
{
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_012(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_013(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_021(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_023(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_031(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_032(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_102(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_103(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_120(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_123(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_130(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_132(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_201(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_203(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_210(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_213(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_230(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_231(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_301(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_302(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_310(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_312(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_320(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_321(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, level, p1FaceStencil);
      
      return;
   } 
}


} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg