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

#include "sor_2D_macroface_P2_update_vertexdofs.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

static void sor_2D_macroface_P2_update_vertexdofs_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * RESTRICT const _data_vertex_stencil_at_vertex, int level, double relax)
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


void sor_2D_macroface_P2_update_vertexdofs(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * RESTRICT const _data_vertex_stencil_at_vertex, int level, double relax)
{
    switch( level )
    {

    default:
        sor_2D_macroface_P2_update_vertexdofs_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edge_stencil_at_vertex, _data_vertexFaceDst, _data_vertexFaceRhs, _data_vertex_stencil_at_vertex, level, relax);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg