
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "communicate_directly_edgedof_cell_to_face_part_1_impl_230.hpp"

namespace hhg {
namespace edgedof {
namespace comm {
namespace generated {

static void communicate_directly_edgedof_cell_to_face_part_1_impl_230_level_any(double const * RESTRICT const _data_edge_cell_src_XYZ, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_X, double * RESTRICT _data_edge_face_dst_gl0_XY, double * RESTRICT _data_edge_face_dst_gl0_XYZ, double * RESTRICT _data_edge_face_dst_gl0_Y, int32_t level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XYZ[(-ctr_1 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) - 1)*(-ctr_1 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 2; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XYZ[(-ctr_1 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) - 1)*(-ctr_1 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 2; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XYZ[(-ctr_1 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) - 1)*(-ctr_1 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 2; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XYZ[(-ctr_1 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) - 1)*(-ctr_1 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 2; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XYZ[(-ctr_1 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) - 1)*(-ctr_1 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 2; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XYZ[(-ctr_1 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) - 1)*(-ctr_1 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 2; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            _data_edge_face_dst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XYZ[(-ctr_1 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 2) - (((-ctr_1 - ctr_2 + (1 << (level)) - 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) - 1)*(-ctr_1 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_1_impl_230(double const * RESTRICT const _data_edge_cell_src_XYZ, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_X, double * RESTRICT _data_edge_face_dst_gl0_XY, double * RESTRICT _data_edge_face_dst_gl0_XYZ, double * RESTRICT _data_edge_face_dst_gl0_Y, int32_t level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_1_impl_230_level_any(_data_edge_cell_src_XYZ, _data_edge_cell_src_Y, _data_edge_cell_src_YZ, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_X, _data_edge_face_dst_gl0_XY, _data_edge_face_dst_gl0_XYZ, _data_edge_face_dst_gl0_Y, level);
        break;
    }
}
    

} // namespace generated
} // namespace comm
} // namespace edgedof
} // namespace hhg