
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "communicate_directly_vertexdof_face_to_cell_impl_103.hpp"

namespace hhg {
namespace vertexdof {
namespace comm {
namespace generated {

static void communicate_directly_vertexdof_face_to_cell_impl_103_level_any(double * RESTRICT _data_p1_cell_dst, double const * RESTRICT const _data_p1_face_src, int32_t level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_cell_dst[-ctr_1 - ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level))] = _data_p1_face_src[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_cell_dst[-ctr_1 - ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level))] = _data_p1_face_src[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)); ctr_1 < (1 << (level)) + 1; ctr_1 += 1)
         {
            _data_p1_cell_dst[-ctr_1 - ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level))] = _data_p1_face_src[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_cell_dst[-ctr_1 - ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level))] = _data_p1_face_src[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_cell_dst[-ctr_1 - ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level))] = _data_p1_face_src[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)); ctr_1 < -ctr_2 + (1 << (level)) + 1; ctr_1 += 1)
         {
            _data_p1_cell_dst[-ctr_1 - ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level))] = _data_p1_face_src[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = (1 << (level)); ctr_2 < (1 << (level)) + 1; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_cell_dst[-ctr_1 - ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level))] = _data_p1_face_src[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}


void communicate_directly_vertexdof_face_to_cell_impl_103(double * RESTRICT _data_p1_cell_dst, double const * RESTRICT const _data_p1_face_src, int32_t level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_face_to_cell_impl_103_level_any(_data_p1_cell_dst, _data_p1_face_src, level);
        break;
    }
}
    

} // namespace generated
} // namespace comm
} // namespace vertexdof
} // namespace hhg