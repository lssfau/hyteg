/*
 * Copyright (c) 2020 Marcus Mohr
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

#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

// using walberla::real_c;
// using walberla::real_t;
// using namespace hyteg;

namespace terraneo {

// type to distinguish different scenarios w.r.t. boundary conditions
typedef enum
{
   BC_NOSLIP,
   BC_FREESLIP,
   BC_MIXED
} bcType;

// ===============
//  SetForceField
// ===============
template < typename vecFuncType >
void setForceField( vecFuncType&                              f,
                    uint_t                                    level,
                    uint_t                                    degree,
                    int                                       order,
                    std::shared_ptr< SphericalHarmonicsTool > sphTool,
                    bcType                                    bmCase )
{
   // safety check
   if ( (uint_t) std::abs( order ) > degree )
   {
      WALBERLA_ABORT( "Spherical harmonics order must be smaller equal to degree!" );
   }

   // first we compute the spherical harmonics contribution which is the same
   // for all three component functions
   std::function< real_t( const Point3D& ) > sphFunc = [sphTool, degree, order]( const Point3D& x ) {
      return sphTool->shconvert_eval( degree, order, x[0], x[1], x[2] );
   };

   typename vecFuncType::VectorComponentType sph( "spherical harmonic", f[0].getStorage(), level, level );
   // feFuncType sph( "spherical harmonic", f[0].getStorage(), level, level );
   sph.interpolate( sphFunc, level, All );

   // set polynomial coefficients for radial component
   real_t                  elfac = real_c( degree * ( degree + 1 ) );
   std::array< real_t, 5 > pc;

   pc[0] = ( 12.0 - elfac ) * ( 2.0 - elfac );
   pc[1] = elfac * ( 6.0 - elfac );
   pc[2] = elfac * ( elfac - 2.0 );
   pc[3] = elfac * ( 2.0 - elfac );
   pc[4] = elfac * ( elfac - 6.0 );

   switch ( bmCase )
   {
   case BC_NOSLIP:
      pc[0] *= 0.25;
      pc[1] *= 1.50;
      pc[2] *= 3.25;
      pc[3] *= 3.00;
      pc[4] *= 1.00;
      break;

   case BC_FREESLIP:
      pc[0] *= 7.0 / 24.0;
      pc[1] *= 1.875;
      pc[2] *= 49.0 / 12.0;
      pc[3] *= 3.5;
      pc[4] *= 1.0;
      break;

   case BC_MIXED:
      pc[0] *= 0.375;
      pc[1] *= 2.125;
      pc[2] *= 4.250;
      pc[3] *= 3.500;
      pc[4] *= 1.000;
      break;
   }

   // x-component
   std::function< real_t( const Point3D& ) > funcX = [pc]( const Point3D& x ) {
      real_t rad = sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      WALBERLA_ASSERT( rad < 2.0 + 1e-12 );
      WALBERLA_ASSERT( rad > 1.0 - 1e-12 );
      real_t iRad  = 1.0 / rad;
      real_t value = ( pc[0] + ( pc[1] + ( pc[2] + ( pc[3] + pc[4] * iRad ) * iRad ) * iRad ) * iRad );

      return value * x[0] * iRad;
   };

   f[0].interpolate( funcX, level, All );
   f[0].multElementwise( {f[0], sph}, level, All );

   // y-component
   std::function< real_t( const Point3D& ) > funcY = [pc]( const Point3D& x ) {
      real_t rad = sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      WALBERLA_ASSERT( rad < 2.0 + 1e-12 );
      WALBERLA_ASSERT( rad > 1.0 - 1e-12 );
      real_t iRad  = 1.0 / rad;
      real_t value = ( pc[0] + ( pc[1] + ( pc[2] + ( pc[3] + pc[4] * iRad ) * iRad ) * iRad ) * iRad );

      return value * x[1] * iRad;
   };

   f[1].interpolate( funcY, level, All );
   f[1].multElementwise( {f[1], sph}, level, All );

   // z-component
   std::function< real_t( const Point3D& ) > funcZ = [pc]( const Point3D& x ) {
      real_t rad = sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
      WALBERLA_ASSERT( rad < 2.0 + 1e-12 );
      WALBERLA_ASSERT( rad > 1.0 - 1e-12 );
      real_t iRad  = 1.0 / rad;
      real_t value = ( pc[0] + ( pc[1] + ( pc[2] + ( pc[3] + pc[4] * iRad ) * iRad ) * iRad ) * iRad );

      return value * x[2] * iRad;
   };

   f[2].interpolate( funcZ, level, All );
   f[2].multElementwise( {f[2], sph}, level, All );
}

// =======================
//  SetAnalyticalSolution
// =======================
template < typename vecFuncType >
void setAnalyticSolution( vecFuncType&                              u,
                          uint_t                                    level,
                          uint_t                                    degree,
                          int                                       order,
                          std::shared_ptr< SphericalHarmonicsTool > sphTool,
                          bcType                                    bmCase )
{
   // safety check
   if ( (uint_t) std::abs( order ) > degree )
   {
      WALBERLA_ABORT( "Spherical harmonics order must be smaller equal to degree!" );
   }

   // set polynomial coefficients for radial component
   real_t                  elfac = real_c( degree * ( degree + 1 ) );
   std::array< real_t, 5 > pCoeff;

   switch ( bmCase )
   {
   case BC_NOSLIP:
      pCoeff[0] = 1.00;
      pCoeff[1] = -3.00;
      pCoeff[2] = 3.25;
      pCoeff[3] = -1.50;
      pCoeff[4] = 0.25;
      break;

   case BC_FREESLIP:
      pCoeff[0] = 1.00;
      pCoeff[1] = -7.0 / 2.0;
      pCoeff[2] = 49.0 / 12.0;
      pCoeff[3] = -15.0 / 8.0;
      pCoeff[4] = 7.0 / 24.0;
      break;

   case BC_MIXED:
      pCoeff[0] = 1.00;
      pCoeff[1] = -7.0 / 2.0;
      pCoeff[2] = 17.0 / 4.0;
      pCoeff[3] = -17.0 / 8.0;
      pCoeff[4] = 3.0 / 8.0;
      break;
   }

   // define function for evaluating the different components
   uint_t                                    component;
   std::function< real_t( const Point3D& ) > vshFunc = [sphTool, degree, order, &component, pCoeff, elfac]( const Point3D& x ) {
      real_t value = real_c( 0 );
      real_t rad   = sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

      // check validity of radial value
      WALBERLA_ASSERT( rad < 2.0 + 1e-12 );
      WALBERLA_ASSERT( rad > 1.0 - 1e-12 );

      real_t poly = ( ( ( pCoeff[4] * rad + pCoeff[3] ) * rad + pCoeff[2] ) * rad + pCoeff[1] ) * rad + pCoeff[0];

      real_t polyDeriv = ( ( 4.0 * pCoeff[4] * rad + 3.0 * pCoeff[3] ) * rad + 2.0 * pCoeff[2] ) * rad + pCoeff[1];

      // Add Y_(l,m)^0 contribution
      value = sphTool->evalVSH( degree, order, x[0], x[1], x[2], 0, component ) * poly * elfac / ( rad * rad );

      // Add Y_(l,m)^1 contribution
      value += sphTool->evalVSH( degree, order, x[0], x[1], x[2], 1, component ) * polyDeriv * sqrt( elfac ) / rad;

      return value;
   };

   // interpolate component functions
   for ( component = 0; component < 3; ++component )
      u[component].interpolate( vshFunc, level, All );
}

} // namespace terraneo
