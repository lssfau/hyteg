
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsEdgeToEdgeMacroFace2D.hpp"

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

static void assign_2D_macroface_edgedof_3_rhs_functions_level_2(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_3(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_4(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_5(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_6(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_7(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_8(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_9(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_10(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_11(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_12(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_13(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_14(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}

static void assign_2D_macroface_edgedof_3_rhs_functions_level_any(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2, int64_t level)
{
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_16 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_17 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_18 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = xi_16 + xi_17 + xi_18;
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
      {
         const double xi_41 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_42 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_43 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_44 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_45 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_46 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = xi_41 + xi_42 + xi_43;
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = xi_44 + xi_45 + xi_46;
      }
      // bottom right vertex
      for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
      {
         const double xi_23 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_24 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_25 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = xi_23 + xi_24 + xi_25;
      }
   }
   for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_65 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_66 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_61 + xi_62 + xi_63;
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = xi_64 + xi_65 + xi_66;
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
      {
         const double xi_3 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_6 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_7 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_8 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_9 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_10 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_11 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = xi_6 + xi_7 + xi_8;
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = xi_10 + xi_11 + xi_9;
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_51 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_52 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_55 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_56 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_51 + xi_52 + xi_53;
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = xi_54 + xi_55 + xi_56;
      }
   }
   for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_30 = c0*_data_edgeFaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_31 = c1*_data_edgeFaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_32 = c2*_data_edgeFaceSrc2[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32;
      }
   }
   {
      
   }
}


void assign_2D_macroface_edgedof_3_rhs_functions(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2, int64_t level)
{
    switch( level )
    {
    case 2:
        assign_2D_macroface_edgedof_3_rhs_functions_level_2(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 3:
        assign_2D_macroface_edgedof_3_rhs_functions_level_3(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 4:
        assign_2D_macroface_edgedof_3_rhs_functions_level_4(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 5:
        assign_2D_macroface_edgedof_3_rhs_functions_level_5(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 6:
        assign_2D_macroface_edgedof_3_rhs_functions_level_6(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 7:
        assign_2D_macroface_edgedof_3_rhs_functions_level_7(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 8:
        assign_2D_macroface_edgedof_3_rhs_functions_level_8(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 9:
        assign_2D_macroface_edgedof_3_rhs_functions_level_9(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 10:
        assign_2D_macroface_edgedof_3_rhs_functions_level_10(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 11:
        assign_2D_macroface_edgedof_3_rhs_functions_level_11(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 12:
        assign_2D_macroface_edgedof_3_rhs_functions_level_12(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 13:
        assign_2D_macroface_edgedof_3_rhs_functions_level_13(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    case 14:
        assign_2D_macroface_edgedof_3_rhs_functions_level_14(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2);
        break;
    default:
        assign_2D_macroface_edgedof_3_rhs_functions_level_any(_data_edgeFaceDst, _data_edgeFaceSrc0, _data_edgeFaceSrc1, _data_edgeFaceSrc2, c0, c1, c2, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg