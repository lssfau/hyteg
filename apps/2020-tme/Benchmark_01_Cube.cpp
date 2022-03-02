/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "Helpers.hpp"

namespace hyteg {
namespace tme_benchmarks {

std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& x ) { return -4 * std::cos( 4 * x[2] ); };
std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& x ) { return 8 * std::cos( 8 * x[0] ); };
std::function< real_t( const hyteg::Point3D& ) > exactW = []( const hyteg::Point3D& x ) { return -2 * std::cos( 2 * x[1] ); };
std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& x ) {
   return std::sin( 4 * x[0] ) * std::sin( 8 * x[1] ) * std::sin( 2 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > rhsU = []( const hyteg::Point3D& x ) {
   return 4 * std::sin( 8 * x[1] ) * std::sin( 2 * x[2] ) * std::cos( 4 * x[0] ) - 64 * std::cos( 4 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > rhsV = []( const hyteg::Point3D& x ) {
   return 8 * std::sin( 4 * x[0] ) * std::sin( 2 * x[2] ) * std::cos( 8 * x[1] ) + 512 * std::cos( 8 * x[0] );
};
std::function< real_t( const hyteg::Point3D& ) > rhsW = []( const hyteg::Point3D& x ) {
   return 2 * std::sin( 4 * x[0] ) * std::sin( 8 * x[1] ) * std::cos( 2 * x[2] ) - 8 * std::cos( 2 * x[1] );
};

void benchmark( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petScManager( &argc, &argv );

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./Benchmark_01_Cube.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const std::string discretizationString = mainConf.getParameter< std::string >( "discretization" );

   const uint_t minLevel        = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel        = mainConf.getParameter< uint_t >( "maxLevel" );
   const uint_t numEdgesPerSide = mainConf.getParameter< uint_t >( "numEdgesPerSide" );

   const bool        vtk    = mainConf.getParameter< bool >( "vtk" );
   const std::string dbFile = mainConf.getParameter< std::string >( "dbFile" );

   MultigridSettings multigridSettings;
   multigridSettings.preSmooth                 = mainConf.getParameter< uint_t >( "preSmooth" );
   multigridSettings.postSmooth                = mainConf.getParameter< uint_t >( "postSmooth" );
   multigridSettings.incSmooth                 = mainConf.getParameter< uint_t >( "incSmooth" );
   multigridSettings.fmgInnerIterations        = mainConf.getParameter< uint_t >( "fmgInnerIterations" );
   multigridSettings.numCycles                 = mainConf.getParameter< uint_t >( "numCycles" );
   multigridSettings.absoluteResidualTolerance = mainConf.getParameter< real_t >( "absoluteResidualTolerance" );

   SmootherSettings smootherSettings;
   smootherSettings.estimateOmega             = mainConf.getParameter< bool >( "estimateOmega" );
   smootherSettings.omega                     = mainConf.getParameter< real_t >( "omega" );
   smootherSettings.omegaEstimationLevel      = mainConf.getParameter< uint_t >( "omegaEstimationLevel" );
   smootherSettings.omegaEstimationIterations = mainConf.getParameter< uint_t >( "omegaEstimationIterations" );
   smootherSettings.numGSVelocity             = mainConf.getParameter< uint_t >( "numGSVelocity" );
   smootherSettings.symmGSVelocity            = mainConf.getParameter< bool >( "symmGSVelocity" );

   CoarseGridSettings coarseGridSettings;
   coarseGridSettings.absoluteResidualTolerance = mainConf.getParameter< real_t >( "coarseGridAbsoluteResidualTolerance" );
   coarseGridSettings.maxIterations             = mainConf.getParameter< uint_t >( "maxIterations" );
   coarseGridSettings.solverType                = mainConf.getParameter< uint_t >( "coarseGridSolverType" );

   const uint_t normCalculationLevelIncrement = mainConf.getParameter< uint_t >( "normCalculationLevelIncrement" );
   const bool   solveWithCoarseGridSolverOnEachFMGLevel =
       mainConf.getParameter< bool >( "solveWithCoarseGridSolverOnEachFMGLevel" );

   Discretization discretization = Discretization::P2_P1;
   if ( discretizationString == "p1p1" )
   {
      discretization = Discretization::P1_P1;
   }

   Point3D leftBottom3D( { 0, 0, 0 } );

   auto meshInfo =
       MeshInfo::meshSymmetricCuboid( leftBottom3D, Point3D( { 1, 1, 1 } ), numEdgesPerSide, numEdgesPerSide, numEdgesPerSide );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   auto onBoundary = []( const Point3D& ) { return true; };
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 1, onBoundary );
   setupStorage->setMeshBoundaryFlagsInner( 0, true );

   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   solve( storage,
          discretization,
          true,
          exactU,
          exactV,
          exactW,
          exactP,
          rhsU,
          rhsV,
          rhsW,
          minLevel,
          maxLevel,
          multigridSettings,
          smootherSettings,
          coarseGridSettings,
          multigridSettings,
          smootherSettings,
          coarseGridSettings,
          true,
          true,
          false,
          normCalculationLevelIncrement,
          solveWithCoarseGridSolverOnEachFMGLevel,
          vtk,
          "Benchmark_01_Cube",
          false,
          dbFile,
          dbFile );
}

} // namespace tme_benchmarks
} // namespace hyteg

int main( int argc, char** argv )
{
   hyteg::tme_benchmarks::benchmark( argc, argv );
}
