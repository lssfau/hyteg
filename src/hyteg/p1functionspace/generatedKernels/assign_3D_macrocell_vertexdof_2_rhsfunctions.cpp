
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "assign_3D_macrocell_vertexdof_2_rhsfunctions.hpp"

namespace hyteg {
namespace vertexdof {
namespace macrocell {
namespace generated {

static void assign_3D_macrocell_vertexdof_2_rhs_functions_level_any(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceSrc0, double * RESTRICT _data_p1FaceSrc1, double c0, double c1, int32_t level)
{
   for (int ctr_3 = 1; ctr_3 < (1 << (level)); ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level)); ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_3 = c0*_data_p1FaceSrc0[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_4 = c1*_data_p1FaceSrc1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            _data_p1FaceDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))] = xi_3 + xi_4;
         }
      }
   }
}


void assign_3D_macrocell_vertexdof_2_rhs_functions(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceSrc0, double * RESTRICT _data_p1FaceSrc1, double c0, double c1, int32_t level)
{
    switch( level )
    {

    default:
        assign_3D_macrocell_vertexdof_2_rhs_functions_level_any(_data_p1FaceDst, _data_p1FaceSrc0, _data_p1FaceSrc1, c0, c1, level);
        break;
    }
}
    

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg