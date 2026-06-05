/*
 * Copyright (c) 2026 Marcus Mohr.
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/forms/form_hyteg_generated/p3/p3_diffusion_affine_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p3/p3_mass_affine_qe.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p3functionspace/P3VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

void printSeparator()
{
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------------------" );
}

void testFormsUnitTriangle()
{
   printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "*** GREY TRIANGLE ***" );
   WALBERLA_LOG_INFO_ON_ROOT( "Testing instantiation of form objects:" );

   WALBERLA_LOG_INFO_ON_ROOT( "-> mass form" );
   forms::p3_mass_affine_qe formMass;

   WALBERLA_LOG_INFO_ON_ROOT( "-> diffusion form" );
   forms::p3_diffusion_affine_q4 formDiffusion;

   printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "Testing evaluation of forms:" );

   // set triangle coordinates
   std::array< Point3D, 3 > coords;
   coords[0] = Point3D{ real_c( 0 ), real_c( 0 ), real_c( 0 ) };
   coords[1] = Point3D{ real_c( 1 ), real_c( 0 ), real_c( 0 ) };
   coords[2] = Point3D{ real_c( 0 ), real_c( 1 ), real_c( 0 ) };

   Matrix< real_t, 10, 10 > elMat;

   Matrix< real_t, 10, 10 > ctrlMatMass;
   ctrlMatMass << .5654761905e-2, .8184523810e-3, .8184523810e-3, .2008928571e-2, .1339285714e-2, .1339285714e-2, .2008928571e-2,
       0., 0., .2678571429e-2, .8184523810e-3, .5654761905e-2, .8184523810e-3, .1339285714e-2, .2008928571e-2, 0., 0.,
       .2008928571e-2, .1339285714e-2, .2678571429e-2, .8184523810e-3, .8184523810e-3, .5654761905e-2, 0., 0., .2008928571e-2,
       .1339285714e-2, .1339285714e-2, .2008928571e-2, .2678571429e-2, .2008928571e-2, .1339285714e-2, 0., .4017857143e-1,
       -.4017857143e-2, -.1004464286e-1, -.1406250000e-1, -.1004464286e-1, .2008928571e-1, .1205357143e-1, .1339285714e-2,
       .2008928571e-2, 0., -.4017857143e-2, .4017857143e-1, .2008928571e-1, -.1004464286e-1, -.1406250000e-1, -.1004464286e-1,
       .1205357143e-1, .1339285714e-2, 0., .2008928571e-2, -.1004464286e-1, .2008928571e-1, .4017857143e-1, -.4017857143e-2,
       -.1004464286e-1, -.1406250000e-1, .1205357143e-1, .2008928571e-2, 0., .1339285714e-2, -.1406250000e-1, -.1004464286e-1,
       -.4017857143e-2, .4017857143e-1, .2008928571e-1, -.1004464286e-1, .1205357143e-1, 0., .2008928571e-2, .1339285714e-2,
       -.1004464286e-1, -.1406250000e-1, -.1004464286e-1, .2008928571e-1, .4017857143e-1, -.4017857143e-2, .1205357143e-1, 0.,
       .1339285714e-2, .2008928571e-2, .2008928571e-1, -.1004464286e-1, -.1406250000e-1, -.1004464286e-1, -.4017857143e-2,
       .4017857143e-1, .1205357143e-1, .2678571429e-2, .2678571429e-2, .2678571429e-2, .1205357143e-1, .1205357143e-1,
       .1205357143e-1, .1205357143e-1, .1205357143e-1, .1205357143e-1, .1446428571;

   formMass.integrateAll( coords, elMat );
   WALBERLA_LOG_INFO_ON_ROOT( "-> mass form\n" << elMat );

   Matrix< real_t, 10, 10 > difference = ( elMat - ctrlMatMass );
   real_t                   diffNorm   = difference.norm();
   WALBERLA_LOG_INFO_ON_ROOT( "-> Frobenius norm of difference to control is " << diffNorm );
   WALBERLA_CHECK_LESS( diffNorm, real_c( 1e-10 ) );

   Matrix< real_t, 10, 10 > ctrlMatStiff;
   ctrlMatStiff << .8499999999, -.8749999995e-1, -.8749999997e-1, -.7499999979e-1, -.6374999998, -.6374999998, -.7499999980e-1,
       .3749999999, .3749999998, -.3403e-9, -.8749999995e-1, .4249999998, 0., .3749999992e-1, -.3750000001e-1, .3374999998,
       .3750000003e-1, -.3750000003e-1, -.6749999997, .940e-10, -.8749999997e-1, 0., .4249999999, .3750000000e-1, .3374999999,
       -.3750000003e-1, .3749999990e-1, -.6749999999, -.3750000000e-1, .770e-10, -.7499999979e-1, .3749999992e-1, .3750000000e-1,
       3.374999999, .3374999997, .3374999994, -.6750000001, .3375000001, -1.687499999, -2.024999999, -.6374999998,
       -.3750000001e-1, .3374999999, .3374999997, 3.375000000, -.2554e-9, .3374999993, -1.687499999, .14002e-9, -2.024999999,
       -.6374999998, .3374999998, -.3750000003e-1, .3374999994, -.2554e-9, 3.374999999, .3374999998, .1206e-9, -1.687499999,
       -2.024999999, -.7499999980e-1, .3750000003e-1, .3749999990e-1, -.6750000001, .3374999993, .3374999998, 3.374999999,
       -1.687499999, .3374999999, -2.024999999, .3749999999, -.3750000003e-1, -.6749999999, .3375000001, -1.687499999, .1206e-9,
       -1.687499999, 3.374999999, .3611e-10, -.5328e-9, .3749999998, -.6749999997, -.3750000000e-1, -1.687499999, .14002e-9,
       -1.687499999, .3374999999, .3611e-10, 3.374999999, -.6091e-9, -.3403e-9, .940e-10, .770e-10, -2.024999999, -2.024999999,
       -2.024999999, -2.024999999, -.5328e-9, -.6091e-9, 8.099999998;

   formDiffusion.integrateAll( coords, elMat );
   WALBERLA_LOG_INFO_ON_ROOT( "-> diffusion form\n" << elMat );

   difference = ( elMat - ctrlMatStiff );
   diffNorm   = difference.norm();
   WALBERLA_LOG_INFO_ON_ROOT( "-> Frobenius norm of difference to control is " << diffNorm );
   WALBERLA_CHECK_LESS( diffNorm, real_c( 1e-8 ) );
}

void testFormsBlueTriangle()
{
   forms::p3_diffusion_affine_q4 formDiffusion;

   printSeparator();
   WALBERLA_LOG_INFO_ON_ROOT( "*** BLUE TRIANGLE ***" );
   WALBERLA_LOG_INFO_ON_ROOT( "Testing evaluation of forms:" );

   // set triangle coordinates
   std::array< Point3D, 3 > coords;
   coords[0] = Point3D{ real_c( 1 ), real_c( 0 ), real_c( 0 ) };
   coords[1] = Point3D{ real_c( 0 ), real_c( 1 ), real_c( 0 ) };
   coords[2] = Point3D{ real_c( 1 ), real_c( 1 ), real_c( 0 ) };

   Matrix< real_t, 10, 10 > elMat;

   Matrix< real_t, 10, 10 > ctrlMatStiff;
   ctrlMatStiff << .4250000000, .1615417985e-18, -.8749999997e-1, -.3749999991e-1, -.6749999999, .3750000010e-1, -.3749999991e-1,
       .3374999999, .3749999991e-1, -.15887e-9, .1615417985e-18, .4249999998, -.8750000000e-1, -.6749999998, -.3750000000e-1,
       .3750000000e-1, .3375000000, -.3750000000e-1, .3750000000e-1, .1536398597e-13, -.8749999997e-1, -.8750000000e-1,
       .8499999999, .3750000000, .3749999999, -.7500000003e-1, -.6375000000, -.6374999999, -.7500000000e-1, .770e-10,
       -.3749999991e-1, -.6749999998, .3750000000, 3.374999999, -.6879e-10, .3374999997, -1.687500000, -.2802093018e-10,
       -1.687499999, .3174e-9, -.6749999999, -.3750000000e-1, .3749999999, -.6879e-10, 3.375000000, -1.687500000, -.2427e-9,
       -1.687500000, .3375000001, .3704e-9, .3750000010e-1, .3750000000e-1, -.7500000003e-1, .3374999997, -1.687500000,
       3.375000000, .3374999999, .3375000001, -.6749999996, -2.025000000, -.3749999991e-1, .3375000000, -.6375000000,
       -1.687500000, -.2427e-9, .3374999999, 3.374999999, .2088e-9, .3375000000, -2.024999999, .3374999999, -.3750000000e-1,
       -.6374999999, -.2802093018e-10, -1.687500000, .3375000001, .2088e-9, 3.375000000, .3375000000, -2.025000000,
       .3749999991e-1, .3750000000e-1, -.7500000000e-1, -1.687499999, .3375000001, -.6749999996, .3375000000, .3375000000,
       3.374999999, -2.025000000, -.15887e-9, .1536398597e-13, .770e-10, .3174e-9, .3704e-9, -2.025000000, -2.024999999,
       -2.025000000, -2.025000000, 8.0999999;

   formDiffusion.integrateAll( coords, elMat );
   WALBERLA_LOG_INFO_ON_ROOT( "-> diffusion form\n" << elMat );

   Matrix< real_t, 10, 10 > difference = ( elMat - ctrlMatStiff );
   real_t                   diffNorm   = difference.norm();

   WALBERLA_LOG_INFO_ON_ROOT( "-> Frobenius norm of difference to control is " << diffNorm );
   WALBERLA_CHECK_LESS( diffNorm, real_c( 2e-7 ) );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::testFormsUnitTriangle();
   hyteg::testFormsBlueTriangle();

   hyteg::printSeparator();

   return EXIT_SUCCESS;
}
