
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "communicate_directly_edgedof_cell_to_face_part_2_impl_302.hpp"

namespace hhg {
namespace edgedof {
namespace comm {
namespace generated {

static void communicate_directly_edgedof_cell_to_face_part_2_impl_302_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int32_t level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_302(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int32_t level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_302_level_any(_data_edge_cell_src_X, _data_edge_cell_src_XY, _data_edge_cell_src_XZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    

} // namespace generated
} // namespace comm
} // namespace edgedof
} // namespace hhg