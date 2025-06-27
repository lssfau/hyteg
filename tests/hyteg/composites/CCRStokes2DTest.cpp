/*
 * Copyright (c) 2017-2025 Marcus Mohr.
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

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"

#include "hyteg/composites/CCRStokesFunction.hpp"
#include "hyteg/composites/CCRStokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/dg1functionspace/DG1Operator.hpp"
#include "hyteg/eigen/EigenSparseDirectSolver.hpp"
#include "hyteg/forms/form_hyteg_dg/DG1MassFormAffine.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "manufactured_solutions/StokesAnalyticalExpressions.hpp"
#include "mixed_operator/P2P1TaylorHoodStokesOperator.hpp"
#include "mixed_operator/VectorMassOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

namespace hyteg {

void objectInstantionTest()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running \"objectInstantionTest\"" );

   uint_t level = 2;

   MeshInfo meshInfo = MeshInfo::meshRectangle(
       Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), hyteg::MeshInfo::CRISSCROSS, 2, 2 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   CCRStokesFunction< real_t > func( "CCR", storage, level, level );
   CCRStokesOperator           oper( storage, level, level );
}

template < typename func_t, typename massOper_t >
void normalisePressureToZeroMean( func_t& pressure, uint_t level )
{
   const std::shared_ptr< PrimitiveStorage > storage = pressure.getStorage();
   massOper_t                                massOperator( storage, level, level );
   func_t                                    aux( "auxilliary function", storage, level, level );
   func_t                                    one( "'ones function'", storage, level, level );
   one.interpolate( real_c( 1 ), level, All );

   massOperator.apply( pressure, aux, level, All, Replace );
   real_t integralValue = aux.dotGlobal( one, level );

   pressure.assign( { real_c( 1 ), -integralValue }, { pressure, one }, level, All );
}

template < typename func_t, typename massOper_t >
real_t computeH0Norm( func_t& feFunction, uint_t level )
{
   const std::shared_ptr< PrimitiveStorage > storage = feFunction.getStorage();
   massOper_t                                massOperator( storage, level, level );

   func_t aux( "auxilliary function", storage, level, level, BoundaryCondition::createAllInnerBC() );

   massOperator.apply( feFunction, aux, level, All, Replace );
   real_t integralValue = aux.dotGlobal( feFunction, level );

   return std::sqrt( integralValue );
}

template < typename stokesOper_t, typename massOperVelocity_t, typename massOperPressure_t >
void solveStokesProblem( uint_t level, bool doVTKOutput, std::string prefix = "" )
{
   // The test uses example D.3 for Steady-State flow problems
   // from John, "Finite Element Methods for Incompressible Flow Problems"
   //             _
   //            |  x^2 (1-x)^4 y^2 (1-y) (3-5y)
   // u = 1000 * |
   //            |_ -2x (1-x)^3 (1-3x) y^3 (1-y)^2
   //
   //
   // p = pi^2 ( x y^3 cos( 2 pi x^2 y ) - x^2 y sin( 2 pi x y ) + 1/8
   //
   // The problem is posed in the unit square with no-slip boundary conditions
   // for velocity.

   using stokesFunc_t   = typename stokesOper_t::srcType;
   using velocityFunc_t = typename massOperVelocity_t::srcType;
   using pressureFunc_t = typename massOperPressure_t::srcType;

   using benchmark_t     = hyteg::manufactured_solutions::Stokes2D::BenchmarkType;
   benchmark_t benchmark = manufactured_solutions::Stokes2D::BenchmarkType::JOHN_D3;
   using expType         = manufactured_solutions::ExpressionType;

   // prepare domain and mesh
   MeshInfo meshInfo = MeshInfo::meshRectangle(
       Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( 1 ), real_c( 1 ) ), hyteg::MeshInfo::CRISSCROSS, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   stokesFunc_t sol_analytic( "Analytic Solution", storage, level, level );

   sol_analytic.uvw().interpolate( { manufactured_solutions::Stokes2D::get( benchmark, expType::VELOCITY_X ),
                                     manufactured_solutions::Stokes2D::get( benchmark, expType::VELOCITY_Y ) },
                                   level );
   sol_analytic.p().interpolate( manufactured_solutions::Stokes2D::get( benchmark, expType::PRESSURE ), level );

   stokesFunc_t rhs_analytic( "Analytic RHS", storage, level, level );
   rhs_analytic.uvw().interpolate( { manufactured_solutions::Stokes2D::get( benchmark, expType::RHS_X ),
                                     manufactured_solutions::Stokes2D::get( benchmark, expType::RHS_Y ) },
                                   level );

   // setup FE problem
   stokesOper_t stokesOperator( storage, level, level );
   stokesFunc_t discrete_rhs( "Discrete RHS", storage, level, level );
   stokesFunc_t discrete_sol( "Discrete Solution", storage, level, level );

   // note: problem has no-slip boundary conditions, so nothing to do for this

   massOperVelocity_t massOperator( storage, level, level );
   massOperator.apply( rhs_analytic.uvw(), discrete_rhs.uvw(), level, All );

#ifdef HYTEG_BUILD_WITH_PETSC
   PETScLUSolver< stokesOper_t > petscLU( storage, level );
   petscLU.solve( stokesOperator, discrete_sol, discrete_rhs, level );
#endif

   // "normalise" pressure to zero mean
   normalisePressureToZeroMean< pressureFunc_t, massOperPressure_t >( discrete_sol.p(), level );

   // compute a simple measure for difference to analytical solutions
   stokesFunc_t difference( "Differences", storage, level, level );
   difference.assign( { real_c( 1 ), real_c( -1 ) }, { sol_analytic, discrete_sol }, level, All );

   real_t vErrorNorm = computeH0Norm< velocityFunc_t, massOperVelocity_t >( difference.uvw(), level );
   real_t pErrorNorm = computeH0Norm< pressureFunc_t, massOperPressure_t >( difference.p(), level );

   WALBERLA_LOG_INFO_ON_ROOT( "vErrorNorm = " << vErrorNorm );
   WALBERLA_LOG_INFO_ON_ROOT( "pErrorNorm = " << pErrorNorm );

   // export data for post-processing
   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", prefix, storage );
      vtkOutput.add( sol_analytic );
      vtkOutput.add( rhs_analytic );
      vtkOutput.add( discrete_sol );
      vtkOutput.add( discrete_rhs );
      vtkOutput.add( difference );
      vtkOutput.write( level );
   }

   WALBERLA_CHECK_LESS_EQUAL( vErrorNorm, real_c( 3e-4 ) );
   WALBERLA_CHECK_LESS_EQUAL( pErrorNorm, real_c( 1.3e-1 ) );
}

void manufacturedSolutionTest( uint_t level, bool doVTKOutput = false, bool runTaylorHoodToo = false )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running \"manufacturedSolutionTest\"" );

   using DG1MassOperator = DG1Operator< DG1MassFormAffine >;
   solveStokesProblem< CCRStokesOperator, P2PlusBubbleElementwiseVectorMassOperator, DG1MassOperator >(
       level, doVTKOutput, "CCRStokes2DTest_with_CCR" );

   // for comparison we might want to solve the same problem with a P2-P1 element
   if ( runTaylorHoodToo )
   {
      solveStokesProblem< P2P1TaylorHoodStokesOperator, P2ConstantVectorMassOperator, P1ConstantMassOperator >(
          level, doVTKOutput, "CCRStokes2DTest_with_TH" );
   }
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   objectInstantionTest();

#ifdef HYTEG_BUILD_WITH_PETSC
   PETScManager petscManager( &argc, &argv );
   manufacturedSolutionTest( 4 );
#endif

   return EXIT_SUCCESS;
}
