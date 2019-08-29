
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "communicate_directly_edgedof_face_to_cell_impl_031.hpp"

namespace hhg {
namespace edgedof {
namespace comm {
namespace generated {

static void communicate_directly_edgedof_face_to_cell_impl_031_level_any(double * RESTRICT _data_edge_cell_dst_X, double * RESTRICT _data_edge_cell_dst_XZ, double * RESTRICT _data_edge_cell_dst_Z, double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, int32_t level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_cell_dst_Z[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_XZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_X[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_cell_dst_Z[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_XZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_X[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_cell_dst_Z[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_XZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_X[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_cell_dst_Z[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_XZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_X[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_cell_dst_Z[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_XZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_X[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_cell_dst_Z[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_XZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_X[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_cell_dst_Z[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_XZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_cell_dst_X[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))] = _data_edge_face_src_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}


void communicate_directly_edgedof_face_to_cell_impl_031(double * RESTRICT _data_edge_cell_dst_X, double * RESTRICT _data_edge_cell_dst_XZ, double * RESTRICT _data_edge_cell_dst_Z, double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, int32_t level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_face_to_cell_impl_031_level_any(_data_edge_cell_dst_X, _data_edge_cell_dst_XZ, _data_edge_cell_dst_Z, _data_edge_face_src_X, _data_edge_face_src_XY, _data_edge_face_src_Y, level);
        break;
    }
}
    

} // namespace generated
} // namespace comm
} // namespace edgedof
} // namespace hhg