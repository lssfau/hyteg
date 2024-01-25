/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#include "core/timing/TimingJSON.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

/*
 * This benchmark meassures the time for several P2 functions on a macro face
 */
int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto              timingTree = std::make_shared< walberla::WcTimingTree >();
   walberla::WcTimer timer;

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if( env.config() == nullptr )
   {
      auto defaultFile = "./P2SolverBenchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   } else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t level    = mainConf.getParameter< uint_t >( "level" );
   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );

   const uint_t cgIterations       = mainConf.getParameter< uint_t >( "cgIterations" );
   const uint_t minresIterations   = mainConf.getParameter< uint_t >( "minresIterations" );
   const uint_t vCycles            = mainConf.getParameter< uint_t >( "vCycles" );
   const uint_t preSmoothingSteps  = mainConf.getParameter< uint_t >( "preSmoothingSteps" );
   const uint_t postSmoothingSteps = mainConf.getParameter< uint_t >( "postSmoothingSteps" );

   const std::string meshFile = mainConf.getParameter< std::string >( "mesh" );

   const bool checkError = mainConf.getParameter< bool >( "checkError" );

   LIKWID_MARKER_THREADINIT;

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   auto storageInfo    = storage->getGlobalInfo();
   auto numVertexDoFs  = numberOfGlobalDoFs< VertexDoFFunctionTag >( *storage, level );
   auto numEdgeDoFs    = numberOfGlobalDoFs< EdgeDoFFunctionTag >( *storage, level );
   auto numP2DoFsTotal = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );

   WALBERLA_LOG_DEVEL_ON_ROOT( "" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "=================================" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "====== P2 Solver Benchmark ======" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "=================================" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );
   WALBERLA_LOG_INFO_ON_ROOT( "mesh:                     " << meshFile );
   WALBERLA_LOG_INFO_ON_ROOT( "refinement level:         " << level );
   WALBERLA_LOG_INFO_ON_ROOT( "# vertexdofs:             " << numVertexDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "# edgedofs:               " << numEdgeDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "# total dofs:             " << numP2DoFsTotal );
   WALBERLA_LOG_INFO_ON_ROOT( "iterations CG:            " << cgIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "iterations MinRes:        " << minresIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "iterations MG (V-cycles): " << vCycles );
   WALBERLA_LOG_INFO_ON_ROOT( "MG smoothing (pre/post):  "
                              << "(" << preSmoothingSteps << ", " << postSmoothingSteps << ")" );

   WALBERLA_LOG_INFO_ON_ROOT( "=== Starting measurements ===" );

   P2Function< real_t > u( "u", storage, minLevel, level );
   P2Function< real_t > f( "f", storage, minLevel, level );
   P2Function< real_t > r( "r", storage, minLevel, level );
   P2Function< real_t > error( "error", storage, minLevel, level );
   P2Function< real_t > exact( "exact", storage, minLevel, level );
   P2Function< real_t > tmp( "tmp", storage, minLevel, level );

   hyteg::P2ConstantLaplaceOperator A( storage, minLevel, level );
   hyteg::P2ConstantMassOperator    M( storage, level, level );

   VTKOutput vtkOutput( ".", "P2SolverBenchmark", storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( r );
   vtkOutput.add( error );
   vtkOutput.add( exact );

   std::function< real_t( const hyteg::Point3D& ) > exactSolution = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > rhsFunction = []( const hyteg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exactSolution, level, DirichletBoundary );
   exact.interpolate( exactSolution, level );
   tmp.interpolate( rhsFunction, level );
   M.apply( tmp, f, level, hyteg::All );

   if( checkError )
   {
      vtkOutput.write( level, 0 );
   }

   ////////
   // CG //
   ////////

   if( cgIterations > 0 )
   {
      u.interpolate( real_c( 0 ), level, Inner );
      CGSolver< P2ConstantLaplaceOperator > cgSolver( storage, level, level, cgIterations );

      LIKWID_MARKER_START( "cg solver" );
      timer.reset();
      cgSolver.solve( A, u, f, level );
      timer.end();
      LIKWID_MARKER_STOP( "cg solver" );
      WALBERLA_LOG_INFO_ON_ROOT( "cg solver:        " << timer.last() );

      if( checkError )
      {
         tmp.interpolate( real_c( 0 ), level );
         A.apply( u, tmp, level, Inner );
         r.assign( {1.0, -1.0}, {f, tmp}, level );
         error.assign( {1.0, -1.0}, {u, exact}, level );
         const real_t residualL2 = std::sqrt( r.dotGlobal( r, level ) / real_c( numP2DoFsTotal ) );
         const real_t errorL2    = std::sqrt( error.dotGlobal( error, level ) / real_c( numP2DoFsTotal ) );
         WALBERLA_LOG_INFO_ON_ROOT( "L2 residual: " << residualL2 );
         WALBERLA_LOG_INFO_ON_ROOT( "L2 error:    " << errorL2 );
         vtkOutput.write( level, 1 );
      }
   }

   ////////////
   // MinRes //
   ////////////

   if( minresIterations > 0 )
   {
      u.interpolate( real_c( 0 ), level, Inner );
      MinResSolver< P2ConstantLaplaceOperator > minresSolver( storage, level, level, minresIterations );

      LIKWID_MARKER_START( "minres solver" );
      timer.reset();
      minresSolver.solve( A, u, f, level );
      timer.end();
      LIKWID_MARKER_STOP( "minres solver" );
      WALBERLA_LOG_INFO_ON_ROOT( "minres solver:    " << timer.last() );

      if( checkError )
      {
         tmp.interpolate( real_c( 0 ), level );
         A.apply( u, tmp, level, Inner );
         r.assign( {1.0, -1.0}, {f, tmp}, level );
         error.assign( {1.0, -1.0}, {u, exact}, level );
         const real_t residualL2 = std::sqrt( r.dotGlobal( r, level ) / real_c( numP2DoFsTotal ) );
         const real_t errorL2    = std::sqrt( error.dotGlobal( error, level ) / real_c( numP2DoFsTotal ) );
         WALBERLA_LOG_INFO_ON_ROOT( "L2 residual: " << residualL2 );
         WALBERLA_LOG_INFO_ON_ROOT( "L2 error:    " << errorL2 );
         vtkOutput.write( level, 2 );
      }
   }

   ////////
   // MG //
   ////////

   if( vCycles > 0 )
   {
      u.interpolate( real_c( 0 ), level, Inner );
      auto smoother             = std::make_shared< GaussSeidelSmoother< P2ConstantLaplaceOperator > >();
      auto coarseGridSolver     = std::make_shared< CGSolver< P2ConstantLaplaceOperator > >( storage, minLevel, minLevel );
      auto restrictionOperator  = std::make_shared< P2toP2QuadraticRestriction >();
      auto prolongationOperator = std::make_shared< P2toP2QuadraticProlongation >();
      GeometricMultigridSolver< P2ConstantLaplaceOperator > gmgSolver( storage,
                                                                       smoother,
                                                                       coarseGridSolver,
                                                                       restrictionOperator,
                                                                       prolongationOperator,
                                                                       minLevel,
                                                                       level,
                                                                       preSmoothingSteps,
                                                                       postSmoothingSteps,
                                                                       0 );

      LIKWID_MARKER_START( "multigrid solver" );
      timer.reset();
      for( uint_t vCycle = 0; vCycle < vCycles; vCycle++ )
      {
         gmgSolver.solve( A, u, f, level );
      }
      timer.end();
      LIKWID_MARKER_STOP( "multigrid solver" );
      WALBERLA_LOG_INFO_ON_ROOT( "multigrid solver: " << timer.last() );

      if( checkError )
      {
         tmp.interpolate( real_c( 0 ), level );
         A.apply( u, tmp, level, Inner );
         r.assign( {1.0, -1.0}, {f, tmp}, level );
         error.assign( {1.0, -1.0}, {u, exact}, level );
         const real_t residualL2 = std::sqrt( r.dotGlobal( r, level ) / real_c( numP2DoFsTotal ) );
         const real_t errorL2    = std::sqrt( error.dotGlobal( error, level ) / real_c( numP2DoFsTotal ) );
         WALBERLA_LOG_INFO_ON_ROOT( "L2 residual: " << residualL2 );
         WALBERLA_LOG_INFO_ON_ROOT( "L2 error:    " << errorL2 );
         vtkOutput.write( level, 3 );
      }
   }

   auto timingTreeReducedWithRemainder = timingTree->getReduced().getCopyWithRemainder();
   WALBERLA_LOG_INFO_ON_ROOT( timingTreeReducedWithRemainder );

   nlohmann::json ttjson = nlohmann::json( timingTreeReducedWithRemainder );
   std::ofstream  o( "P2SolverBenchmark.json" );
   o << ttjson;
   o.close();

   LIKWID_MARKER_CLOSE;
}
