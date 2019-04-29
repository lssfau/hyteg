
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsVertexToVertexMacroFace3D.hpp"

namespace hhg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void apply_3D_macroface_vertexdof_to_vertexdof_replace_130_130(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, double const * RESTRICT const _data_p1FaceSrc_gl1, int64_t level, std::map< walberla::uint_t, std::map< hhg::indexing::IndexIncrement, double > > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[0][{ 0, 1, -1 }];
   const double xi_2 = p1FaceStencil[0][{ -1, 1, -1 }];
   const double xi_3 = p1FaceStencil[0][{ 0, 1, 0 }];
   const double xi_4 = p1FaceStencil[0][{ -1, 1, 0 }];
   const double xi_5 = p1FaceStencil[1][{ 0, 1, -1 }];
   const double xi_6 = p1FaceStencil[1][{ -1, 1, -1 }];
   const double xi_7 = p1FaceStencil[1][{ 0, 1, 0 }];
   const double xi_8 = p1FaceStencil[1][{ -1, 1, 0 }];
   const double xi_9 = p1FaceStencil[0][{ 1, 0, -1 }];
   const double xi_10 = p1FaceStencil[1][{ 1, 0, -1 }];
   const double xi_11 = p1FaceStencil[0][{ 0, 0, -1 }];
   const double xi_12 = p1FaceStencil[1][{ 0, 0, -1 }];
   const double xi_13 = p1FaceStencil[0][{ 1, 0, 0 }];
   const double xi_14 = p1FaceStencil[1][{ 1, 0, 0 }];
   const double xi_15 = p1FaceStencil[0][{ 0, 0, 0 }];
   const double xi_16 = p1FaceStencil[1][{ 0, 0, 0 }];
   const double xi_17 = p1FaceStencil[0][{ -1, 0, 0 }];
   const double xi_18 = p1FaceStencil[1][{ -1, 0, 0 }];
   const double xi_19 = p1FaceStencil[0][{ 0, 0, 1 }];
   const double xi_20 = p1FaceStencil[1][{ 0, 0, 1 }];
   const double xi_21 = p1FaceStencil[0][{ -1, 0, 1 }];
   const double xi_22 = p1FaceStencil[1][{ -1, 0, 1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_25 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_40 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_41 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = xi_5*_data_p1FaceSrc_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_43 = xi_6*_data_p1FaceSrc_gl1[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_44 = xi_7*_data_p1FaceSrc_gl1[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_45 = xi_8*_data_p1FaceSrc_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_9*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = xi_10*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_27 = xi_11*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_28 = xi_12*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_29 = xi_13*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_30 = xi_14*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_31 = xi_15*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_16*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_33 = xi_17*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_34 = xi_18*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_35 = xi_19*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_37 = xi_20*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_38 = xi_21*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_39 = xi_22*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46;
      }
   }
}

static void apply_3D_macroface_vertexdof_to_vertexdof_replace_130_132(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, double const * RESTRICT const _data_p1FaceSrc_gl1, int64_t level, std::map< walberla::uint_t, std::map< hhg::indexing::IndexIncrement, double > > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[0][{ 0, 1, -1 }];
   const double xi_2 = p1FaceStencil[0][{ -1, 1, -1 }];
   const double xi_3 = p1FaceStencil[0][{ 0, 1, 0 }];
   const double xi_4 = p1FaceStencil[0][{ -1, 1, 0 }];
   const double xi_5 = p1FaceStencil[1][{ 0, 0, -1 }];
   const double xi_6 = p1FaceStencil[1][{ -1, 1, -1 }];
   const double xi_7 = p1FaceStencil[1][{ 0, -1, 0 }];
   const double xi_8 = p1FaceStencil[1][{ -1, 0, 0 }];
   const double xi_9 = p1FaceStencil[0][{ 1, 0, -1 }];
   const double xi_10 = p1FaceStencil[1][{ 1, 0, -1 }];
   const double xi_11 = p1FaceStencil[0][{ 0, 0, -1 }];
   const double xi_12 = p1FaceStencil[1][{ 0, 1, -1 }];
   const double xi_13 = p1FaceStencil[0][{ 1, 0, 0 }];
   const double xi_14 = p1FaceStencil[1][{ 1, -1, 0 }];
   const double xi_15 = p1FaceStencil[0][{ 0, 0, 0 }];
   const double xi_16 = p1FaceStencil[1][{ 0, 0, 0 }];
   const double xi_17 = p1FaceStencil[0][{ -1, 0, 0 }];
   const double xi_18 = p1FaceStencil[1][{ -1, 1, 0 }];
   const double xi_19 = p1FaceStencil[0][{ 0, 0, 1 }];
   const double xi_20 = p1FaceStencil[1][{ 0, -1, 1 }];
   const double xi_21 = p1FaceStencil[0][{ -1, 0, 1 }];
   const double xi_22 = p1FaceStencil[1][{ -1, 0, 1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_25 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_40 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_41 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = xi_5*_data_p1FaceSrc_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_43 = xi_6*_data_p1FaceSrc_gl1[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_44 = xi_7*_data_p1FaceSrc_gl1[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_45 = xi_8*_data_p1FaceSrc_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_9*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = xi_10*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_27 = xi_11*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_28 = xi_12*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_29 = xi_13*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_30 = xi_14*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_31 = xi_15*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_16*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_33 = xi_17*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_34 = xi_18*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_35 = xi_19*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_37 = xi_20*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_38 = xi_21*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_39 = xi_22*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46;
      }
   }
}

static void apply_3D_macroface_vertexdof_to_vertexdof_replace_132_130(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, double const * RESTRICT const _data_p1FaceSrc_gl1, int64_t level, std::map< walberla::uint_t, std::map< hhg::indexing::IndexIncrement, double > > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[0][{ 0, 0, -1 }];
   const double xi_2 = p1FaceStencil[0][{ -1, 1, -1 }];
   const double xi_3 = p1FaceStencil[0][{ 0, -1, 0 }];
   const double xi_4 = p1FaceStencil[0][{ -1, 0, 0 }];
   const double xi_5 = p1FaceStencil[1][{ 0, 1, -1 }];
   const double xi_6 = p1FaceStencil[1][{ -1, 1, -1 }];
   const double xi_7 = p1FaceStencil[1][{ 0, 1, 0 }];
   const double xi_8 = p1FaceStencil[1][{ -1, 1, 0 }];
   const double xi_9 = p1FaceStencil[0][{ 1, 0, -1 }];
   const double xi_10 = p1FaceStencil[1][{ 1, 0, -1 }];
   const double xi_11 = p1FaceStencil[0][{ 0, 1, -1 }];
   const double xi_12 = p1FaceStencil[1][{ 0, 0, -1 }];
   const double xi_13 = p1FaceStencil[0][{ 1, -1, 0 }];
   const double xi_14 = p1FaceStencil[1][{ 1, 0, 0 }];
   const double xi_15 = p1FaceStencil[0][{ 0, 0, 0 }];
   const double xi_16 = p1FaceStencil[1][{ 0, 0, 0 }];
   const double xi_17 = p1FaceStencil[0][{ -1, 1, 0 }];
   const double xi_18 = p1FaceStencil[1][{ -1, 0, 0 }];
   const double xi_19 = p1FaceStencil[0][{ 0, -1, 1 }];
   const double xi_20 = p1FaceStencil[1][{ 0, 0, 1 }];
   const double xi_21 = p1FaceStencil[0][{ -1, 0, 1 }];
   const double xi_22 = p1FaceStencil[1][{ -1, 0, 1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_25 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_40 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_41 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = xi_5*_data_p1FaceSrc_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_43 = xi_6*_data_p1FaceSrc_gl1[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_44 = xi_7*_data_p1FaceSrc_gl1[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_45 = xi_8*_data_p1FaceSrc_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_9*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = xi_10*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_27 = xi_11*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_28 = xi_12*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_29 = xi_13*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_30 = xi_14*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_31 = xi_15*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_16*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_33 = xi_17*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_34 = xi_18*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_35 = xi_19*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_37 = xi_20*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_38 = xi_21*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_39 = xi_22*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46;
      }
   }
}

static void apply_3D_macroface_vertexdof_to_vertexdof_replace_132_132(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, double const * RESTRICT const _data_p1FaceSrc_gl1, int64_t level, std::map< walberla::uint_t, std::map< hhg::indexing::IndexIncrement, double > > p1FaceStencil)
{
   const double xi_1 = p1FaceStencil[0][{ 0, 0, -1 }];
   const double xi_2 = p1FaceStencil[0][{ -1, 1, -1 }];
   const double xi_3 = p1FaceStencil[0][{ 0, -1, 0 }];
   const double xi_4 = p1FaceStencil[0][{ -1, 0, 0 }];
   const double xi_5 = p1FaceStencil[1][{ 0, 0, -1 }];
   const double xi_6 = p1FaceStencil[1][{ -1, 1, -1 }];
   const double xi_7 = p1FaceStencil[1][{ 0, -1, 0 }];
   const double xi_8 = p1FaceStencil[1][{ -1, 0, 0 }];
   const double xi_9 = p1FaceStencil[0][{ 1, 0, -1 }];
   const double xi_10 = p1FaceStencil[1][{ 1, 0, -1 }];
   const double xi_11 = p1FaceStencil[0][{ 0, 1, -1 }];
   const double xi_12 = p1FaceStencil[1][{ 0, 1, -1 }];
   const double xi_13 = p1FaceStencil[0][{ 1, -1, 0 }];
   const double xi_14 = p1FaceStencil[1][{ 1, -1, 0 }];
   const double xi_15 = p1FaceStencil[0][{ 0, 0, 0 }];
   const double xi_16 = p1FaceStencil[1][{ 0, 0, 0 }];
   const double xi_17 = p1FaceStencil[0][{ -1, 1, 0 }];
   const double xi_18 = p1FaceStencil[1][{ -1, 1, 0 }];
   const double xi_19 = p1FaceStencil[0][{ 0, -1, 1 }];
   const double xi_20 = p1FaceStencil[1][{ 0, -1, 1 }];
   const double xi_21 = p1FaceStencil[0][{ -1, 0, 1 }];
   const double xi_22 = p1FaceStencil[1][{ -1, 0, 1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_25 = xi_1*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_36 = xi_2*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_40 = xi_3*_data_p1FaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_41 = xi_4*_data_p1FaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = xi_5*_data_p1FaceSrc_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_43 = xi_6*_data_p1FaceSrc_gl1[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_44 = xi_7*_data_p1FaceSrc_gl1[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_45 = xi_8*_data_p1FaceSrc_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = xi_9*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_26 = xi_10*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_27 = xi_11*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_28 = xi_12*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_29 = xi_13*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_30 = xi_14*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_31 = xi_15*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = xi_16*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_33 = xi_17*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_34 = xi_18*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_35 = xi_19*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_37 = xi_20*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_38 = xi_21*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_39 = xi_22*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46;
      }
   }
}

void apply_3D_macroface_vertexdof_to_vertexdof_replace(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, double const * RESTRICT const _data_p1FaceSrc_gl1, int64_t level, int64_t neighbor_cell_0_local_vertex_id_0, int64_t neighbor_cell_0_local_vertex_id_1, int64_t neighbor_cell_0_local_vertex_id_2, int64_t neighbor_cell_1_local_vertex_id_0, int64_t neighbor_cell_1_local_vertex_id_1, int64_t neighbor_cell_1_local_vertex_id_2, std::map< walberla::uint_t, std::map< hhg::indexing::IndexIncrement, double > > p1FaceStencil)
{
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      apply_3D_macroface_vertexdof_to_vertexdof_replace_130_130(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, _data_p1FaceSrc_gl1, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      apply_3D_macroface_vertexdof_to_vertexdof_replace_130_132(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, _data_p1FaceSrc_gl1, level, p1FaceStencil);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      apply_3D_macroface_vertexdof_to_vertexdof_replace_132_130(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, _data_p1FaceSrc_gl1, level, p1FaceStencil);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      apply_3D_macroface_vertexdof_to_vertexdof_replace_132_132(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceSrc_gl0, _data_p1FaceSrc_gl1, level, p1FaceStencil);
      
      return;
   } 
}


} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg