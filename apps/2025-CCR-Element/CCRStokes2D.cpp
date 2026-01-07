/*
 * Copyright (c) 2025 Marcus Mohr.
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

#include "hyteg/ccrfunctionspace/CCRStokesFunction.hpp"
#include "hyteg/ccrfunctionspace/CCRStokesOperator.hpp"
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

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This app only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using benchmark_t = hyteg::manufactured_solutions::Stokes2D::BenchmarkType;

namespace hyteg {

struct ResultsTabulator
{
   void insertTestResult( uint_t level, std::pair< real_t, real_t > errorNorms, bool resultForTH = false )
   {
      if ( !resultForTH )
      {
         testResultsCCR_[level] = errorNorms;
      }
      else
      {
         testResultsTH_[level] = errorNorms;
      }
   }

   void outputTable( uint_t minLevel, uint_t maxLevel )
   {
      printTable( minLevel, maxLevel, "   Results for the Conforming Crouzeix Raviart Element   ", testResultsCCR_ );

      if ( testResultsTH_.size() > 0 )
      {
         printTable( minLevel, maxLevel, "        Results for the P2-P1 Taylor-Hood Element        ", testResultsTH_ );
      }
   }

   void printTable( uint_t                                           minLevel,
                    uint_t                                           maxLevel,
                    std::string                                      headerText,
                    std::map< uint_t, std::pair< real_t, real_t > >& results )
   {
      std::stringstream table;

      table << "-------------------------------------------------------------\n";
      table << "| " << headerText << " |\n";
      table << "|-----------------------------------------------------------|\n";
      table << "| level | velocity error | factor | pressure error | factor |\n";
      table << "|-----------------------------------------------------------|\n";

      real_t vErrorOld = real_c( 0 );
      real_t pErrorOld = real_c( 0 );

      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         real_t vErrorCur = results[level].first;
         real_t pErrorCur = results[level].second;

         // level info
         table << "| " << std::setw( 3 ) << level << "   |  ";

         // velocity info
         table << std::scientific << std::setprecision( 6 ) << vErrorCur << "  | ";
         if ( level == minLevel )
         {
            table << " ---- ";
         }
         else
         {
            table << std::setprecision( 3 ) << std::setw( 6 ) << std::fixed << vErrorOld / vErrorCur;
         }
         vErrorOld = vErrorCur;

         // pressure info
         table << " |  " << std::scientific << std::setprecision( 6 ) << pErrorCur << "  | ";
         if ( level == minLevel )
         {
            table << " ---- ";
         }
         else
         {
            table << std::setprecision( 4 ) << std::fixed << pErrorOld / pErrorCur;
         }
         pErrorOld = pErrorCur;

         table << " |" << std::endl;
      }

      table << "-------------------------------------------------------------\n";

      WALBERLA_LOG_INFO_ON_ROOT( "" << table.str() );
   }

   std::map< uint_t, std::pair< real_t, real_t > > testResultsCCR_;
   std::map< uint_t, std::pair< real_t, real_t > > testResultsTH_;
};

std::shared_ptr< PrimitiveStorage > genStorageForUnitSquare( uint_t numSubIntervals )
{
   MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( real_c( 0 ), real_c( 0 ) ),
                                                Point2D( real_c( 1 ), real_c( 1 ) ),
                                                hyteg::MeshInfo::CRISSCROSS,
                                                numSubIntervals,
                                                numSubIntervals );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   return std::make_shared< PrimitiveStorage >( setupStorage, 0 );
}

template < typename func_t, typename massOper_t >
void normalisePressureToZeroMean( func_t& pressure, uint_t level )
{
   const std::shared_ptr< PrimitiveStorage > storage = pressure.getStorage();
   massOper_t                                massOperator( storage, level, level );

   func_t aux( "auxilliary function", storage, level, level, BoundaryCondition::createAllInnerBC() );
   func_t one( "'ones function'", storage, level, level, BoundaryCondition::createAllInnerBC() );

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
auto solveStokesProblem( benchmark_t                          benchmark,
                         std::shared_ptr< PrimitiveStorage >& storage,
                         uint_t                               level,
                         bool                                 doVTKOutput,
                         std::string                          prefix = "" )
{
   using stokesFunc_t   = typename stokesOper_t::srcType;
   using pressureFunc_t = typename massOperPressure_t::srcType;
   using velocityFunc_t = typename massOperVelocity_t::srcType;

   stokesFunc_t sol_analytic( "Analytic Solution", storage, level, level );

   using expType = manufactured_solutions::ExpressionType;

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

   // EigenSparseDirectSolver< stokesOper_t > EigenLU( storage, level );
   // EigenLU.solve( stokesOperator, discrete_sol, discrete_rhs, level );

   PETScLUSolver< stokesOper_t > petscLU( storage, level );
   petscLU.solve( stokesOperator, discrete_sol, discrete_rhs, level );

   // "normalise" pressure to zero mean
   normalisePressureToZeroMean< pressureFunc_t, massOperPressure_t >( discrete_sol.p(), level );

   // compute difference to analytical solutions
   stokesFunc_t difference( "Differences", storage, level, level );
   difference.assign( { real_c( 1 ), real_c( -1 ) }, { sol_analytic, discrete_sol }, level, All );

   // NOTE: We currently compute the error norms on the computational level; should improve this

   real_t vErrorNorm = computeH0Norm< velocityFunc_t, massOperVelocity_t >( difference.uvw(), level );
   real_t pErrorNorm = computeH0Norm< pressureFunc_t, massOperPressure_t >( difference.p(), level );

   WALBERLA_LOG_INFO_ON_ROOT( "  -> " << prefix << ": velocity error = " << std::scientific << vErrorNorm );
   WALBERLA_LOG_INFO_ON_ROOT( "  -> " << prefix << ": pressure error = " << std::scientific << pErrorNorm );

   // export data for post-processing
   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( "output", prefix, storage );
      vtkOutput.add( sol_analytic );
      vtkOutput.add( rhs_analytic );
      vtkOutput.add( discrete_sol );
      vtkOutput.add( discrete_rhs.uvw() );
      vtkOutput.add( difference );
      vtkOutput.write( level );
   }

   return std::make_pair( vErrorNorm, pErrorNorm );
}

std::shared_ptr< ResultsTabulator > runBenchmark( benchmark_t                          benchmark,
                                                  std::shared_ptr< PrimitiveStorage >& storage,
                                                  uint_t                               minLevel,
                                                  uint_t                               maxLevel,
                                                  bool                                 compareTaylorHood,
                                                  bool                                 doVTKOutput )
{
   using DG1MassOperator = DG1Operator< DG1MassFormAffine >;

   std::pair< real_t, real_t > errorNorms;

   auto benchmarkResults = std::make_shared< ResultsTabulator >();

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "* level = " << level );

      errorNorms = solveStokesProblem< CCRStokesOperator, P2PlusBubbleElementwiseVectorMassOperator, DG1MassOperator >(
          benchmark, storage, level, doVTKOutput, "CCRStokes2DTest_with_CCR" );
      benchmarkResults->insertTestResult( level, errorNorms );

      if ( compareTaylorHood )
      {
         errorNorms = solveStokesProblem< P2P1TaylorHoodStokesOperator, P2ConstantVectorMassOperator, P1ConstantMassOperator >(
             benchmark, storage, level, doVTKOutput, "CCRStokes2DTest_with_TH" );
         benchmarkResults->insertTestResult( level, errorNorms, true );
      }
   }

   return benchmarkResults;
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   // Setup enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );
   walberla::MPIManager::instance()->useWorldComm();

   // Read steering parameters (if no config file was given on command line load default one)
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( walberlaEnv.config() == nullptr )
   {
      auto defaultFile = "./CCRStokes2D.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = walberlaEnv.config();
   }
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   // Fire up PETSc to use as linear solver
   PETScManager petscManager( &argc, &argv );

   // Select benchmark case according to parameter file
   benchmark_t benchmark;
   std::string bmCase = parameters.getParameter< std::string >( "benchmark" );
   if ( bmCase.compare( "DONEA_HUERTA" ) == 0 )
   {
      benchmark = manufactured_solutions::Stokes2D::BenchmarkType::DONEA_HUERTA;
   }
   else if ( bmCase.compare( "JOHN_D3" ) == 0 )
   {
      benchmark = manufactured_solutions::Stokes2D::BenchmarkType::JOHN_D3;
   }
   else
   {
      WALBERLA_ABORT( "Benchmark case '" << bmCase << "' not supported" );
   }

   // Extract further parameters
   uint_t minLevel          = parameters.getParameter< uint_t >( "minLevel" );
   uint_t maxLevel          = parameters.getParameter< uint_t >( "maxLevel" );
   uint_t numSubIntervals   = parameters.getParameter< uint_t >( "numSubIntervals" );
   bool   doVTKOutput       = parameters.getParameter< bool >( "outputVTK" );
   bool   compareTaylorHood = parameters.getParameter< bool >( "compareTaylorHood" );

   std::shared_ptr< PrimitiveStorage > storage = genStorageForUnitSquare( numSubIntervals );

   WALBERLA_LOG_INFO_ON_ROOT( "Running Benchmark: " << bmCase );
   std::shared_ptr< ResultsTabulator > benchm =
       runBenchmark( benchmark, storage, minLevel, maxLevel, compareTaylorHood, doVTKOutput );

   benchm->outputTable( minLevel, maxLevel );

   return EXIT_SUCCESS;
}
