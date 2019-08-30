
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "communicate_directly_vertexdof_cell_to_face_impl_023.hpp"

namespace hyteg {
namespace vertexdof {
namespace comm {
namespace generated {

static void communicate_directly_vertexdof_cell_to_face_impl_023_level_any(double const * RESTRICT const _data_p1_cell_src, double * RESTRICT _data_p1_face_dst_gl0, int32_t level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
   }
}


void communicate_directly_vertexdof_cell_to_face_impl_023(double const * RESTRICT const _data_p1_cell_src, double * RESTRICT _data_p1_face_dst_gl0, int32_t level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_023_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}
    

} // namespace generated
} // namespace comm
} // namespace vertexdof
} // namespace hyteg