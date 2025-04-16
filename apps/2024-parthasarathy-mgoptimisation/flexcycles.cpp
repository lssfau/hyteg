/*
 * Copyright (c) 2024 Dinesh Parthasarathy.
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

#include "core/Environment.h"
#include "core/Hostname.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/BuildInfo.hpp"
#include "hyteg/Git.hpp"
#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP2Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP1Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVersion.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/EmptySolver.hpp"
#include "hyteg/solvers/FlexibleMultigridSolver.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/SORSmoother.hpp"
#include "hyteg/solvers/SymmetricGaussSeidelSmoother.hpp"
#include "hyteg/solvers/SymmetricSORSmoother.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg_operators/operators/div_k_grad/P1ElementwiseDivKGrad.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGrad.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "constant_stencil_operator/P2ConstantOperator.hpp"

namespace hyteg {

using walberla::int64_c;
using walberla::int8_t;
using walberla::int_c;
using walberla::math::pi;

// solution and rhs: 2d poisson
std::function< real_t( const hyteg::Point3D& ) > exact;
std::function< real_t( const hyteg::Point3D& ) > rhs;

std::function< real_t( const hyteg::Point3D& ) > exactConstk = []( const hyteg::Point3D& x ) {
   return cos( pi * x[0] ) - sin( 2.0 * pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > rhsConstk = []( const hyteg::Point3D& x ) {
   return pi * pi * cos( pi * x[0] ) - 4.0 * pi * pi * sin( 2.0 * pi * x[1] );
};

enum class MeshType
{
   SQUARE,
   CUBE,
   SYMMETRIC_CUBE,
   SPHERICAL_SHELL,
   T_DOMAIN,
   SNAKE
};

const std::map< std::string, MeshType > meshTypeStrings = {
    { "square", MeshType::SQUARE },
    { "cube", MeshType::CUBE },
    { "symmetricCube", MeshType::SYMMETRIC_CUBE },
    { "sphericalShell", MeshType::SPHERICAL_SHELL },
    { "tDomain", MeshType::T_DOMAIN },
    { "snake", MeshType::SNAKE },
};

// calculate error and residual
template < typename Function, typename ElementWiseOperator, typename MassOperator >
void calculateErrorAndResidual( const uint_t&              level,
                                const ElementWiseOperator& A,
                                const MassOperator&,
                                const Function& u,
                                const Function& f,
                                const Function& uExact,
                                const Function& error,
                                const Function& residual,
                                const Function& tmp,
                                long double&    LInfError,
                                long double&    L2Error,
                                long double&    LInfResidual,
                                long double&    L2Residual )
{
   error.assign( { 1.0, -1.0 }, { uExact, u }, level, All );

   tmp.interpolate( real_c( 0 ), level, All );
   A.apply( u, tmp, level, Inner );
   residual.assign( { 1.0, -1.0 }, { f, tmp }, level, All );

   auto num = numberOfGlobalDoFs< typename Function::Tag >( *u.getStorage(), level );

   LInfError    = error.getMaxDoFMagnitude( level );
   LInfResidual = residual.getMaxDoFMagnitude( level );

   L2Error    = std::sqrt( error.dotGlobal( error, level, Inner ) / (long double) ( num ) );
   L2Residual = std::sqrt( residual.dotGlobal( residual, level, Inner ) / (long double) ( num ) );
}

// calculate discretization error
template < typename Function, typename ElementWiseOperator, typename MassOperator >
void calculateDiscretizationError( const std::shared_ptr< PrimitiveStorage >& storage,
                                   const ElementWiseOperator&                 A,
                                   const uint_t&                              level,
                                   long double&                               l2DiscretizationError )
{
   Function u( "u", storage, level, level );
   Function f( "f", storage, level, level );

   Function uExact( "uExact", storage, level, level );
   Function residual( "residual", storage, level, level );
   Function error( "error", storage, level, level );
   Function tmp( "tmp", storage, level, level );

   MassOperator M( storage, level, level );

   u.interpolate( exact, level, DirichletBoundary );
   uExact.interpolate( exact, level, All );

   tmp.interpolate( rhs, level, All );
   M.apply( tmp, f, level, All );

#ifdef HYTEG_BUILD_WITH_PETSC
   auto solver = std::make_shared< PETScLUSolver< ElementWiseOperator > >( storage, level );
#else
   auto solver = std::make_shared< CGSolver< ElementWiseOperator > >( storage, level, level );
#endif
   solver->solve( A, u, f, level );

   long double L2Error;
   long double l2Residual;
   long double L2Residual;
   calculateErrorAndResidual(
       level, A, M, u, f, uExact, error, residual, tmp, L2Error, l2DiscretizationError, l2Residual, L2Residual );
}

template < typename Function,
           typename ElementWiseOperator,
           typename MassOperator,
           typename Restriction,
           typename Prolongation,
           typename FMGProlongation >
void MultigridSolve( const std::shared_ptr< PrimitiveStorage >& storage,
                     const ElementWiseOperator&                 A,
                     const uint_t&                              minLevel,
                     const uint_t&                              maxLevel,
                     const uint_t&                              numCycles,
                     const CycleType                            cycleType,
                     const uint_t&                              fmgInnerCycles,
                     const bool&                                flexCycle,
                     const std::vector< int8_t >                cycleStructure,
                     const std::vector< int >                   smootherTypes,
                     const std::vector< real_t >                smootherWeights,
                     const real_t&                              L2residualTolerance,
                     const uint_t&                              preSmoothingSteps,
                     const uint_t&                              postSmoothingSteps,
                     const uint_t&                              coarseGridMaxIterations,
                     const real_t&                              coarseGridResidualTolerance,
                     const bool&                                outputVTK,
                     const bool&                                calcDiscretizationError )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Solving Poisson problem using multigrid..." );
   Function u( "u", storage, minLevel, maxLevel );
   Function f( "f", storage, minLevel, maxLevel );
   Function x0( "x0", storage, minLevel, maxLevel );

   Function uExact( "uExact", storage, minLevel, maxLevel );
   Function residual( "residual", storage, minLevel, maxLevel );
   Function error( "error", storage, minLevel, maxLevel );
   Function tmp( "tmp", storage, minLevel, maxLevel );

   MassOperator M( storage, minLevel, maxLevel );

   //////////////////////////////////////////////
   // Initialize functions and right-hand side //
   //////////////////////////////////////////////

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      u.interpolate( exact, level, DirichletBoundary );
      uExact.interpolate( exact, level, All );

      tmp.interpolate( rhs, level, All );
      M.apply( tmp, f, level, All );
   }

   long double LInfError;
   long double L2Error;
   long double LInfResidual;
   long double L2Residual;
   long double L2ResidualInitial;

   ////////////////////
   // Initialize VTK //
   ////////////////////

   VTKOutput vtkOutput( "vtk", "P2MultigridLaplace", storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( uExact );
   vtkOutput.add( residual );
   vtkOutput.add( error );

   /////////////////////////
   // Misc setup and info //
   /////////////////////////

   long double discretizationError = 0.0;
   if ( calcDiscretizationError )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "l2 discretization error per level:" );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         calculateDiscretizationError< Function, ElementWiseOperator, MassOperator >( storage, A, level, discretizationError );
         WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::scientific << discretizationError );
      }
      WALBERLA_LOG_DEVEL_ON_ROOT( "" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Number of unknowns (including boundary):" )
   uint_t totalDoFs = 0;
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      const uint_t dofsThisLevel = numberOfGlobalDoFs< typename Function::Tag >( *storage, level );
      WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::setw( 15 ) << dofsThisLevel );
      totalDoFs += dofsThisLevel;
   }
   WALBERLA_LOG_INFO_ON_ROOT( " ----------------------------- " );
   WALBERLA_LOG_INFO_ON_ROOT( "  total:    " << std::setw( 15 ) << totalDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   walberla::WcTimer timer;
   double            timeError;
   double            timeVTK;
   double            timeCycle, timeMGSolve = 0.0;

   timer.reset();
   calculateErrorAndResidual( maxLevel, A, M, u, f, uExact, error, residual, tmp, LInfError, L2Error, LInfResidual, L2Residual );
   timer.end();
   timeError         = timer.last();
   L2ResidualInitial = L2Residual;
   timer.reset();
   if ( outputVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }
   timer.end();
   timeVTK = timer.last();

   WALBERLA_LOG_INFO_ON_ROOT(
       " After cycle... ||   LInf error |     L2 error | L2 error reduction ||  LInf res.   |  L2 residual | L2 residual reduction || time cycle [s] | time error calculation [s] | time VTK [s] |" );
   WALBERLA_LOG_INFO_ON_ROOT(
       " ---------------++--------------+--------------+--------------------++--------------+--------------+-----------------------++----------------+----------------------------+--------------|" );
   WALBERLA_LOG_INFO_ON_ROOT( "        initial || " << std::scientific << LInfError << " | " << L2Error << " | "
                                                    << "               --- || " << LInfResidual << " | " << L2Residual
                                                    << " |                   --- ||            --- | " << std::fixed
                                                    << std::setprecision( 2 ) << std::setw( 26 ) << timeError << " | "
                                                    << std::setw( 12 ) << timeVTK << " |" );

   long double avgL2ErrorConvergenceRate      = 0;
   long double avgL2ResidualConvergenceRate   = 0;
   long double avgLInfErrorConvergenceRate    = 0;
   long double avgLInfResidualConvergenceRate = 0;

   long double L2ErrorReduction      = 0;
   long double L2ResidualReduction   = 0;
   long double LInfErrorReduction    = 0;
   long double LInfResidualReduction = 0;

   ///////////
   // Solve //
   ///////////
   auto smoother = std::make_shared< GaussSeidelSmoother< ElementWiseOperator > >();
   std::vector< std::shared_ptr< Solver< ElementWiseOperator > > > smootherList;
   if ( flexCycle )
   {
      uint_t n_nodes = cycleStructure.size() + 1;
      for ( uint_t i = 0; i < n_nodes; i++ )
      {
         switch ( smootherTypes[i] )
         {
         case evostencils::SolverType::SymmetricSOR:
            smootherList.push_back( std::make_shared< SymmetricSORSmoother< ElementWiseOperator > >( smootherWeights[i] ) );
            break;
         case evostencils::SolverType::WeightedJacobi:
            smootherList.push_back( std::make_shared< WeightedJacobiSmoother< ElementWiseOperator > >(
                storage, minLevel, maxLevel, smootherWeights[i] ) );
            break;
         case evostencils::SolverType::GaussSeidel:
            smootherList.push_back( std::make_shared< GaussSeidelSmoother< ElementWiseOperator > >() );
            break;
         case evostencils::SolverType::SymmetricGaussSeidel:
            smootherList.push_back( std::make_shared< SymmetricGaussSeidelSmoother< ElementWiseOperator > >() );
            break;
         case evostencils::SolverType::SOR:
            smootherList.push_back( std::make_shared< SORSmoother< ElementWiseOperator > >( smootherWeights[i] ) );
            break;
         case evostencils::SolverType::Empty:
            smootherList.push_back( std::make_shared< EmptySolver< ElementWiseOperator > >() );
            break;
         }
      }
   }
#ifdef HYTEG_BUILD_WITH_PETSC
   auto coarseGridSolver = std::make_shared< PETScCGSolver< ElementWiseOperator > >(
       storage, minLevel, coarseGridResidualTolerance, 1e-16, coarseGridMaxIterations );
#else
   auto coarseGridSolver = std::make_shared< CGSolver< ElementWiseOperator > >(
       storage, minLevel, minLevel, coarseGridMaxIterations, coarseGridResidualTolerance );
#endif

   auto prolongationOperator = std::make_shared< Prolongation >();
   auto restrictionOperator  = std::make_shared< Restriction >();
   auto multigridflexSolver  = std::make_shared< FlexibleMultigridSolver< ElementWiseOperator > >(
       storage, smootherList, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, cycleStructure );
   auto multigridSolver = std::make_shared< GeometricMultigridSolver< ElementWiseOperator > >( storage,
                                                                                               smoother,
                                                                                               coarseGridSolver,
                                                                                               restrictionOperator,
                                                                                               prolongationOperator,
                                                                                               minLevel,
                                                                                               maxLevel,
                                                                                               preSmoothingSteps,
                                                                                               postSmoothingSteps,
                                                                                               0,
                                                                                               cycleType );

   auto fmgProlongation = std::make_shared< FMGProlongation >();

   auto postCycle = [&]( uint_t currentLevel ) {
      long double _l2Error, _L2Error, _l2Residual, _L2Residual;
      calculateErrorAndResidual(
          currentLevel, A, M, u, f, uExact, error, residual, tmp, _l2Error, _L2Error, _l2Residual, _L2Residual );
      WALBERLA_LOG_INFO_ON_ROOT( "    fmg level " << currentLevel << ": l2 error: " << std::scientific << _l2Error );
   };

   FullMultigridSolver< ElementWiseOperator > fullMultigridSolver(
       storage, multigridSolver, fmgProlongation, minLevel, maxLevel, fmgInnerCycles, postCycle );
   FullMultigridSolver< ElementWiseOperator, FlexibleMultigridSolver< ElementWiseOperator > > fullflexMultigridSolver(
       storage, multigridflexSolver, fmgProlongation, minLevel, maxLevel, fmgInnerCycles, postCycle );

   uint_t numExecutedCycles = 0;
   if ( flexCycle )
      WALBERLA_LOG_INFO_ON_ROOT( "using flexible multigrid solver" );
   for ( uint_t cycle = 1; cycle <= numCycles; cycle++ )
   {
      const long double lastL2Error      = L2Error;
      const long double lastL2Residual   = L2Residual;
      const long double lastLInfError    = LInfError;
      const long double lastLInfResidual = LInfResidual;

      timer.reset();
      if ( cycle == 1 && fmgInnerCycles > 0 )
      {
         if ( flexCycle )
         {
            fullflexMultigridSolver.solve( A, u, f, maxLevel );
         }
         else
         {
            fullMultigridSolver.solve( A, u, f, maxLevel );
         }
      }
      else
      {
         if ( flexCycle )
         {
            multigridflexSolver->solve( A, u, f, maxLevel );
         }
         else
         {
            multigridSolver->solve( A, u, f, maxLevel );
         }
      }
      timer.end();
      timeCycle = timer.last();
      timeMGSolve += timeCycle;

      numExecutedCycles++;

      timer.reset();
      calculateErrorAndResidual(
          maxLevel, A, M, u, f, uExact, error, residual, tmp, LInfError, L2Error, LInfResidual, L2Residual );
      timer.end();
      timeError = timer.last();

      timer.reset();
      if ( outputVTK )
      {
         vtkOutput.write( maxLevel, cycle );
      }
      timer.end();
      timeVTK = timer.last();

      L2ErrorReduction      = L2Error / lastL2Error;
      L2ResidualReduction   = L2Residual / lastL2Residual;
      LInfErrorReduction    = LInfError / lastLInfError;
      LInfResidualReduction = LInfResidual / lastLInfResidual;

      WALBERLA_LOG_INFO_ON_ROOT(
          std::setw( 15 ) << cycle << " || " << std::scientific << LInfError << " | " << L2Error << " | "
                          << "      " << L2ErrorReduction << " || " << LInfResidual << " | " << L2Residual << " |          "
                          << L2ResidualReduction << " || " << std::fixed << std::setprecision( 2 ) << std::setw( 14 ) << timeCycle
                          << " | " << std::setw( 26 ) << timeError << " | " << std::setw( 12 ) << timeVTK << " | "
                          << " | ratio discr.err: " << ( calcDiscretizationError ? LInfError / discretizationError : 0.0 ) );

      avgL2ErrorConvergenceRate += L2ErrorReduction;
      avgL2ResidualConvergenceRate += L2ResidualReduction;
      avgLInfErrorConvergenceRate += LInfErrorReduction;
      avgLInfResidualConvergenceRate += LInfResidualReduction;

      if ( L2Residual < L2residualTolerance )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "L2 residual dropped below tolerance." )
         break;
      }
   }

   avgL2ErrorConvergenceRate /= real_c( numExecutedCycles );
   avgL2ResidualConvergenceRate = std::pow( L2Residual / L2ResidualInitial, 1.0 / real_c( numExecutedCycles ) );

   avgLInfErrorConvergenceRate /= real_c( numExecutedCycles );
   avgLInfResidualConvergenceRate /= real_c( numExecutedCycles );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Average convergence rates:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 error:      " << std::scientific << avgL2ErrorConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual:   " << std::scientific << avgL2ResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - LInf error:    " << std::scientific << avgLInfErrorConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - LInf residual: " << std::scientific << avgLInfResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "Solve Time: " << timeMGSolve );
   WALBERLA_LOG_INFO_ON_ROOT( "Convergence Factor: " << avgL2ResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of Iterations: " << numExecutedCycles );
}

void setup( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "///////////////////////" );
   WALBERLA_LOG_INFO_ON_ROOT( "// Multigrid Flexible Cycles//" );
   WALBERLA_LOG_INFO_ON_ROOT( "///////////////////////" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   buildinfo::printGitInfo();
   buildinfo::printBuildInfo();
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./flexcycles.prm";
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   ////////////////
   // Parameters //
   ////////////////

   // problem parameters
   const uint_t      numProcesses                 = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
   const uint_t      numEdgesPerSide              = mainConf.getParameter< uint_t >( "numEdgesPerSide" );
   const std::string discretization               = mainConf.getParameter< std::string >( "discretization" );
   const real_t      L2residualTolerance          = mainConf.getParameter< real_t >( "L2residualTolerance" );
   const bool        calculateDiscretizationError = mainConf.getParameter< bool >( "calculateDiscretizationError" );
   // multigrid parameters
   const uint_t      numCycles          = mainConf.getParameter< uint_t >( "numCycles" );
   const std::string cycleTypeString    = mainConf.getParameter< std::string >( "cycleType" );
   const uint_t      fmgInnerCycles     = mainConf.getParameter< uint_t >( "fmgInnerCycles" );
   const uint_t      preSmoothingSteps  = mainConf.getParameter< uint_t >( "preSmoothingSteps" );
   const uint_t      postSmoothingSteps = mainConf.getParameter< uint_t >( "postSmoothingSteps" );
   const uint_t      smoothingIncrement = mainConf.getParameter< uint_t >( "smoothingIncrement" );
   uint_t            minLevel           = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel           = mainConf.getParameter< uint_t >( "maxLevel" );
   // coarse grid solver parameters
   const uint_t coarseGridMaxIterations     = mainConf.getParameter< uint_t >( "coarseGridMaxIterations" );
   real_t       coarseGridResidualTolerance = mainConf.getParameter< real_t >( "coarseGridResidualTolerance" );
   const bool   blockLowRank                = mainConf.getParameter< bool >( "blockLowRank" );
   const real_t blockLowRankTolerance       = mainConf.getParameter< real_t >( "blockLowRankTolerance" );
   // flexible cycle parameters
   const bool  flexCycle             = mainConf.getParameter< bool >( "flexCycle" );
   std::string cycleStructureString  = mainConf.getParameter< std::string >( "cycleStructure" );
   std::string smootherTypesString   = mainConf.getParameter< std::string >( "smootherTypes" );
   std::string smootherWeightsString = mainConf.getParameter< std::string >( "smootherWeights" );
   // output parameters
   const std::string outputBaseDirectory  = mainConf.getParameter< std::string >( "outputBaseDirectory" );
   const bool        outputVTK            = mainConf.getParameter< bool >( "outputVTK" );
   const bool        outputTiming         = mainConf.getParameter< bool >( "outputTiming" );
   const bool        outputTimingJSON     = mainConf.getParameter< bool >( "outputTimingJSON" );
   const std::string outputTimingJSONFile = mainConf.getParameter< std::string >( "outputTimingJSONFile" );

   for ( int i = 1; i < argc; ++i )
   {
      if ( strcmp( argv[i], "-cycleStructure" ) == 0 )
      {
         cycleStructureString = argv[i + 1];
         i++;
      }
      else if ( strcmp( argv[i], "-smootherTypes" ) == 0 )
      {
         smootherTypesString = argv[i + 1];
         i++;
      }
      else if ( strcmp( argv[i], "-smootherWeights" ) == 0 )
      {
         smootherWeightsString = argv[i + 1];
         i++;
      }
      else if ( strcmp( argv[i], "-coarseGridResidualTolerance" ) == 0 )
      {
         coarseGridResidualTolerance = atof( argv[i + 1] );
         i++;
      }
      else if ( strcmp( argv[i], "-minLevel" ) == 0 )
      {
         minLevel = static_cast< uint_t >( atoi( argv[i + 1] ) );
         i++;
      }
   }

   const CycleType       cycleType = ( cycleTypeString == "V" ? CycleType::VCYCLE : CycleType::WCYCLE );
   std::vector< int8_t > cycleStructure;
   // parse cycle structure string into vector of int8_t
   if ( cycleStructureString != "" )
   {
      std::stringstream ss( cycleStructureString );
      int               i;
      while ( ss >> i )
      {
         if ( ss.peek() == ',' )
         {
            ss.ignore(); // Skip the comma
         }
         cycleStructure.push_back( static_cast< int8_t >( i ) );
      }
   }

   std::vector< int > smootherTypes;
   // parse smoother types string into vector of int8_t
   if ( smootherTypesString != "" )
   {
      std::stringstream ss( smootherTypesString );
      int               i;
      while ( ss >> i )
      {
         smootherTypes.push_back( i );
         if ( ss.peek() == ',' )
         {
            ss.ignore();
         }
      }
   }

   std::vector< real_t > smootherWeights;
   // parse smoother weights string into vector of real_t
   if ( smootherWeightsString != "" )
   {
      std::stringstream ss( smootherWeightsString );
      real_t            i;
      while ( ss >> i )
      {
         smootherWeights.push_back( i );
         if ( ss.peek() == ',' )
         {
            ss.ignore();
         }
      }
   }

#ifdef HYTEG_BUILD_WITH_PETSC
   PETScManager petscManager( &argc, &argv );
   printPETScVersionNumberString();
   WALBERLA_LOG_INFO_ON_ROOT( "" );
#endif

   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num processes:                           " << numProcesses );
   WALBERLA_LOG_INFO_ON_ROOT( "  - discretization:                          " << discretization );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num edges per side:                      " << numEdgesPerSide );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num cycles:                              " << numCycles );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual tolerance:                   " << L2residualTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "  - min / max level:                         " << minLevel << " / " << maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "  - flexible cycle:                          " << ( flexCycle ? "yes" : "no" ) );
   if ( flexCycle )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - cycle structure:                         " << cycleStructureString );
      WALBERLA_LOG_INFO_ON_ROOT( "  - smoother types:                          " << smootherTypesString );
      WALBERLA_LOG_INFO_ON_ROOT( "  - smoother weights:                        " << smootherWeightsString );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  - cycle type:                              " << cycleTypeString );
      WALBERLA_LOG_INFO_ON_ROOT( "  - pre- / post- / incr-smoothing:           "
                                 << preSmoothingSteps << " / " << postSmoothingSteps << " / " << smoothingIncrement );
   }
   WALBERLA_LOG_INFO_ON_ROOT(
       "  - full multigrid:                          "
       << ( fmgInnerCycles == 0 ? "no" : "yes, inner cycles per level: " + std::to_string( fmgInnerCycles ) ) );

   WALBERLA_LOG_INFO_ON_ROOT( "  - coarse grid max iterations: " << coarseGridMaxIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "  - coarse grid residual tol: " << coarseGridResidualTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "  - BLR:                                     " << ( blockLowRank ? "enabled" : "disabled" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - BLR tolerance:                           " << blockLowRankTolerance );

   WALBERLA_LOG_INFO_ON_ROOT( "  - output base directory:                   " << outputBaseDirectory );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output VTK:                              " << ( outputVTK ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing:                           " << ( outputTiming ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing JSON:                      " << ( outputTimingJSON ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing JSON file:                 " << outputTimingJSONFile );

   exact = exactConstk;
   rhs   = rhsConstk;
   // setup the mesh: 2d square domain with criss-cross mesh
   std::shared_ptr< PrimitiveStorage > storage;
   Point2D                             leftBottom( 0, 0 );
   auto meshInfo = MeshInfo::meshRectangle( leftBottom, Point2D( 1, 1 ), MeshInfo::CRISSCROSS, numEdgesPerSide, numEdgesPerSide );
   SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   storage         = std::make_shared< PrimitiveStorage >( setupStorage );
   auto globalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( globalInfo );
   if ( outputVTK )
   {
      writeDomainPartitioningVTK( storage, "vtk", outputBaseDirectory + "/Domain" );
   }

   if ( discretization == "P1" )
   {
      P1ConstantLaplaceOperator A( storage, minLevel, maxLevel );
      MultigridSolve< P1Function< real_t >,
                      P1ConstantLaplaceOperator,
                      P1ConstantMassOperator,
                      P1toP1LinearRestriction<>,
                      P1toP1LinearProlongation<>,
                      P1toP1QuadraticProlongation >( storage,
                                                     A,
                                                     minLevel,
                                                     maxLevel,
                                                     numCycles,
                                                     cycleType,
                                                     fmgInnerCycles,
                                                     flexCycle,
                                                     cycleStructure,
                                                     smootherTypes,
                                                     smootherWeights,
                                                     L2residualTolerance,
                                                     preSmoothingSteps,
                                                     postSmoothingSteps,
                                                     coarseGridMaxIterations,
                                                     coarseGridResidualTolerance,
                                                     outputVTK,
                                                     calculateDiscretizationError );
   }
   else if ( discretization == "P2" )
   {
      P2ConstantLaplaceOperator A( storage, minLevel, maxLevel );
      MultigridSolve< P2Function< real_t >,
                      P2ConstantLaplaceOperator,
                      P2ConstantMassOperator,
                      P2toP2QuadraticRestriction,
                      P2toP2QuadraticProlongation,
                      P2toP2QuadraticProlongation >( storage,
                                                     A,
                                                     minLevel,
                                                     maxLevel,
                                                     numCycles,
                                                     cycleType,
                                                     fmgInnerCycles,
                                                     flexCycle,
                                                     cycleStructure,
                                                     smootherTypes,
                                                     smootherWeights,
                                                     L2residualTolerance,
                                                     preSmoothingSteps,
                                                     postSmoothingSteps,
                                                     coarseGridMaxIterations,
                                                     coarseGridResidualTolerance,
                                                     outputVTK,
                                                     calculateDiscretizationError );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Discretization " << discretization << " not supported." );
   }

   if ( outputTiming )
   {
      printTimingTree( *storage->getTimingTree() );
   }

   if ( outputTimingJSON )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing JSON timing to " << outputTimingJSONFile )
      writeTimingTreeJSON( *storage->getTimingTree(), outputTimingJSONFile );
      WALBERLA_LOG_INFO_ON_ROOT( "Done writing JSON timing." )
   }
}
} // namespace hyteg

int main( int argc, char** argv )
{
   hyteg::setup( argc, argv );
}
