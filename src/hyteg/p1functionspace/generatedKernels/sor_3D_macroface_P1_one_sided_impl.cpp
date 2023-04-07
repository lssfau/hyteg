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

#include "sor_3D_macroface_P1_one_sided_impl.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void sor_3D_macroface_P1_one_sided_impl_012_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_11*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_8*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_9*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_012(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_012_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_013_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_11*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_5*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_6*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_013(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_013_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_021_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_5*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_8*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_9*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_021(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_021_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_023_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_4*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_7*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_023(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_023_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_031_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_5*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_031(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_031_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_032_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_4*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_6*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_032(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_032_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_102_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_11*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_3*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_9*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_102(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_102_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_103_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_5*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_103(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_103_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_120_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_3*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_9*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_120(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_120_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_123_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
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
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_5*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_123(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_123_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_130_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_3*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_130(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_130_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_132_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
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
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_132(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_132_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_201_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_3*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_5*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_7*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_201(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_201_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_203_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_4*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_6*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_203(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_203_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_210_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_3*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_7*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_210(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_210_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_213_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
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
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_2*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_5*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_213(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_213_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_230_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_11*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_3*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_4*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_6*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_230(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_230_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_231_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
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
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_2*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_231(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_231_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_301_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_5*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_7*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_301(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_301_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_302_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_5*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_6*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_302(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_302_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_310_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_3*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_6*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_310(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_310_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_312_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
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
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_11*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_2*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_312(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_312_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_320_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_13 = 1 / (xi_1);
   const double xi_2 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_3 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_4 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_5 = v2v_cell_stencil_fused_face_0[{ 0, 0, 1 }];
   const double xi_6 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_7 = v2v_cell_stencil_fused_face_0[{ 0, 1, 0 }];
   const double xi_8 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_9 = v2v_cell_stencil_fused_face_0[{ 1, -1, 1 }];
   const double xi_10 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_11 = v2v_cell_stencil_fused_face_0[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_2*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_3*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_4*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_320(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_320_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    
static void sor_3D_macroface_P1_one_sided_impl_321_level_any(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
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
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_13*(-xi_10*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_11*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_2*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_3*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_5*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_6*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_7*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_8*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_9*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P1_one_sided_impl_321(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int level, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P1_one_sided_impl_321_level_any(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg