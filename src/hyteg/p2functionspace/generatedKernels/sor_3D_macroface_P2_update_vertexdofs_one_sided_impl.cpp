/*
 * Copyright (c) 2019 Nils Kohl, Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "sor_3D_macroface_P2_update_vertexdofs_one_sided_impl.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_012_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_13*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_14*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_15*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_17*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_18*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_19*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_24*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_25*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_27*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_28*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_29*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_31*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_33*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_34*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_35*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_38*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_39*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_40*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_42*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_012(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_012_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_013_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_13*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_14*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_15*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_16*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_19*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_23*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_24*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_25*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_27*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_28*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_31*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_32*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_33*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_36*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_37*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_41*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_42*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_013(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_013_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_021_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_12*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_13*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_14*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_15*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_16*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_17*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_21*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_24*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_25*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_28*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_29*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_30*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_34*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_35*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_36*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_38*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_39*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_5*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_edgeFaceDst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_021(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_021_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_023_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_13*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_14*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_15*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_16*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_19*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_24*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_25*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_28*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_29*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_31*_data_edgeFaceDst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_33*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_35*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_37*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_38*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_42*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_5*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_023(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_023_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_031_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_13*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_14*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_15*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_16*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_18*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_24*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_25*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_27*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_28*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_30*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_31*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_36*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_37*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_38*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_41*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_5*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_031(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_031_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_032_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_13*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_14*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_15*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_16*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_18*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_22*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_23*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_24*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_25*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_27*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_28*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_30*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_35*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_36*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_42*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_5*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_6*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_7*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_032(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_032_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_102_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_13*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_17*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_18*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_20*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_21*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_24*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_28*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_31*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_34*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_38*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_39*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_42*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_6*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_8*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_102(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_102_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_103_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_12*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_13*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_16*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_17*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_19*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_20*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_21*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_23*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_24*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_28*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_29*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_36*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_37*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_41*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_42*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_6*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_7*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_103(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_103_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_120_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_13*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_14*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_16*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_17*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_21*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_23*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_24*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_28*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_30*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_34*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_38*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_39*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_5*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_8*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_120(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_120_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_123_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_12*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_13*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_14*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_15*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_16*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_17*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_18*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_21*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_24*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_25*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_27*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_28*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_29*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_30*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_31*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_33*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_36*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_41*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_5*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_123(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_123_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_130_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_13*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_14*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_16*_data_edgeFaceDst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_18*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_24*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_27*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_28*_data_edgeFaceDst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_30*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_31*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_34*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_36*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_37*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_41*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_5*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_7*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_130(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_130_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_132_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_13*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_14*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_15*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_17*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_18*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_19*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_20*_data_edgeFaceDst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_24*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_25*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_27*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_28*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_41*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_42*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_6*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_7*_data_edgeFaceDst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_132(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_132_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_201_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_12*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_13*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_15*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_18*_data_edgeFaceDst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_20*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_21*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_24*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_25*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_28*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_29*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_31*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_34*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_35*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_36*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_39*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_6*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_edgeFaceDst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_201(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_201_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_203_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_13*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_14*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_15*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_17*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_18*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_19*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_20*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_21*_data_edgeFaceDst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_24*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_25*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_28*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_29*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_35*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_37*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_42*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_203(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_203_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_210_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_13*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_14*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_15*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_18*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_19*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_23*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_24*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_25*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_27*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_28*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_29*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_31*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_34*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_35*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_39*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_40*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_5*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_210(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_210_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_213_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_13*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_14*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_16*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_19*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_23*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_24*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_25*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_26*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_28*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_29*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_33*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_36*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_40*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_41*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_5*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_213(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_213_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_230_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_12*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_13*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_14*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_15*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_16*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_17*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_18*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_22*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_24*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_25*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_27*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_28*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_30*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_34*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_35*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_42*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_5*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_7*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_230(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_230_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_231_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_13*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_14*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_17*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_19*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_20*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_24*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_25*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_27*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_28*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_41*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_6*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_7*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_231(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_231_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_301_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_12*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_13*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_14*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_15*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_17*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_19*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_20*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_21*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_24*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_25*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_28*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_29*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_36*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_37*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_38*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_41*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_6*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_8*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_301(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_301_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_302_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_13*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_14*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_17*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_18*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_19*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_20*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_21*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_23*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_24*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_28*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_29*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_35*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_36*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_37*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_42*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_6*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_8*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_302(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_302_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_310_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_13*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_14*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_15*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_16*_data_edgeFaceDst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_19*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_24*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_25*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_27*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_28*_data_edgeFaceDst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_29*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_31*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_32*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_34*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_36*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_37*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_41*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_5*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_8*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_310(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_310_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_312_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_13*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_14*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_15*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_16*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_18*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_19*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_23*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_24*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_26*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_28*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_29*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_31*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_32*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_33*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_34*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_40*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_41*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_42*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_312(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_312_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_320_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_12*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_13*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_14*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_16*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_17*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_18*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_19*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_24*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_26*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_28*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_29*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_30*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_31*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_33*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_34*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_35*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_37*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_40*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_41*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_42*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_5*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_8*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_320(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_320_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_321_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_44 = 1 / (xi_1);
   const double xi_2 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_33 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_34 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_35 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_36 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_37 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_38 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_39 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_40 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_41 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_42 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_44*(-xi_10*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_12*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_13*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_14*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_15*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_16*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_17*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_18*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_19*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_21*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_22*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_23*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_24*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_26*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_27*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_28*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_29*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_30*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_31*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_32*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_33*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_34*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_35*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_36*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_38*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_39*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_40*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_41*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_42*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_5*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_6*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_7*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_321(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil_fused_face_0, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_one_sided_impl_321_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg