/*
 * Copyright (c) 2017-2020 Nils Kohl, Dominik Thoennes.
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
namespace scaling_workshop {

void benchmark( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petScManager( &argc, &argv );

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./Scaling_Workshop_02_Spherical_Shell.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   std::string benchmarkName = "Benchmark_02_Spherical_Shell";

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const std::string discretizationString = mainConf.getParameter< std::string >( "discretization" );
   const uint_t      minLevel             = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel             = mainConf.getParameter< uint_t >( "maxLevel" );
   const uint_t      ntan                 = mainConf.getParameter< uint_t >( "ntan" );
   const uint_t      nrad                 = mainConf.getParameter< uint_t >( "nrad" );
   const real_t      rmin                 = mainConf.getParameter< real_t >( "rmin" );
   const real_t      rmax                 = mainConf.getParameter< real_t >( "rmax" );

   const bool        vtk            = mainConf.getParameter< bool >( "vtk" );
   const std::string dbFile         = mainConf.getParameter< std::string >( "dbFile" );
   const std::string timingFile     = mainConf.getParameter< std::string >( "timingFile" );
   const bool        domainInfoOnly = mainConf.getParameter< bool >( "domainInfoOnly" );

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

   const uint_t scenario = mainConf.getParameter< uint_t >( "scenario" );

   Discretization discretization = Discretization::P2_P1;
   if ( discretizationString == "p1p1" )
   {
      discretization = Discretization::P1_P1;
   }

   std::shared_ptr< PrimitiveStorage > storage;
   {
      auto meshInfo = MeshInfo::meshSphericalShell( ntan, nrad, rmin, rmax, MeshInfo::SHELLMESH_ON_THE_FLY );

      SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      // new code ...
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      // ... replaces
      // auto onBoundary = []( const Point3D& ) { return true; };
      // meshInfo.setMeshBoundaryFlagsByVertexLocation( 1, onBoundary );
      // setupStorage.setMeshBoundaryFlagsInner( 0, true );

      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }

   if ( domainInfoOnly )
   {
      domainInfo( storage, discretization, minLevel, maxLevel, vtk, benchmarkName );
      return;
   }

   std::function< real_t( const hyteg::Point3D& ) > solutionU = []( const hyteg::Point3D& ) -> real_t { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > solutionV = []( const hyteg::Point3D& ) -> real_t { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > solutionW = []( const hyteg::Point3D& ) -> real_t { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > solutionP = []( const hyteg::Point3D& ) -> real_t { return 0; };

   std::function< real_t( const hyteg::Point3D& ) > initialU = []( const hyteg::Point3D& ) -> real_t { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > initialV = []( const hyteg::Point3D& ) -> real_t { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > initialW = []( const hyteg::Point3D& ) -> real_t { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > initialP = []( const hyteg::Point3D& ) -> real_t { return 0; };

   std::function< real_t( const hyteg::Point3D& ) > rhsU = []( const hyteg::Point3D& ) -> real_t { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > rhsV = []( const hyteg::Point3D& ) -> real_t { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > rhsW = []( const hyteg::Point3D& ) -> real_t { return 0; };

   bool projectPressure                 = true;
   bool projectPressureAfterRestriction = true;
   bool RHSisZero                       = false;

   FixedSizeSQLDB db( dbFile );
   db.setConstantEntry( "ntan", ntan );
   db.setConstantEntry( "nrad", nrad );
   db.setConstantEntry( "rmin", rmin );
   db.setConstantEntry( "rmax", rmax );

   WALBERLA_LOG_INFO_ON_ROOT( "########################" )
   WALBERLA_LOG_INFO_ON_ROOT( "### Scaling Workshop ###" )
   WALBERLA_LOG_INFO_ON_ROOT( "########################" )
   WALBERLA_LOG_INFO_ON_ROOT( "# domain: spherical shell" )
   WALBERLA_LOG_INFO_ON_ROOT( "# - ntan: " << ntan )
   WALBERLA_LOG_INFO_ON_ROOT( "# - nrad: " << nrad )
   WALBERLA_LOG_INFO_ON_ROOT( "# - rmin: " << rmin )
   WALBERLA_LOG_INFO_ON_ROOT( "# - rmax: " << rmax )
   if ( scenario == 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "# - scenario 0: u, p = 0, f = 0, initial guess: rand(0, 1) " );

      initialU = []( const hyteg::Point3D& ) { return walberla::math::realRandom(); };
      initialV = []( const hyteg::Point3D& ) { return walberla::math::realRandom(); };
      initialW = []( const hyteg::Point3D& ) { return walberla::math::realRandom(); };
      initialP = []( const hyteg::Point3D& ) { return walberla::math::realRandom(); };

      RHSisZero = true;
   }
   else if ( scenario == 1 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "# - scenario 1: Stokeslet, f = (1, 1, 1), initial guess 0 " );

      solutionU = []( const hyteg::Point3D& p ) -> real_t {
        const Point3D fDirection( {1, 1, 1} );
        const real_t factorU = 1.0 / (8.0 * walberla::math::pi );

        const uint_t i = 0;
         const auto r    = p.norm();
         const auto rInv = 1.0 / r;

         Point3D q;

         for ( uint_t j = 0; j < 3; j++ )
         {
            q[j] = p[i] * p[j] * rInv * rInv * rInv;
         }
         q[i] += rInv;

         const auto result = factorU * q.dot( fDirection );
         return result;
      };

      solutionV = []( const hyteg::Point3D& p ) -> real_t {
        const Point3D fDirection( {1, 1, 1} );
        const real_t factorU = 1.0 / (8.0 * walberla::math::pi );

        const uint_t i = 1;
        const auto r    = p.norm();
        const auto rInv = 1.0 / r;

        Point3D q;

        for ( uint_t j = 0; j < 3; j++ )
        {
           q[j] = p[i] * p[j] * rInv * rInv * rInv;
        }
        q[i] += rInv;

        const auto result = factorU * q.dot( fDirection );
        return result;
      };

      solutionW = []( const hyteg::Point3D& p ) -> real_t {
        const Point3D fDirection( {1, 1, 1} );
        const real_t factorU = 1.0 / (8.0 * walberla::math::pi );

        const uint_t i = 2;
        const auto r    = p.norm();
        const auto rInv = 1.0 / r;

        Point3D q;

        for ( uint_t j = 0; j < 3; j++ )
        {
           q[j] = p[i] * p[j] * rInv * rInv * rInv;
        }
        q[i] += rInv;

        const auto result = factorU * q.dot( fDirection );
        return result;
      };

      solutionP = []( const hyteg::Point3D& p ) -> real_t {
        const Point3D fDirection( {1, 1, 1} );
        const real_t factorP = 1.0 / (4.0 * walberla::math::pi );

        const auto r    = p.norm();
        const auto rInv = 1.0 / r;

        Point3D q = p;
        q *= rInv * rInv * rInv;

        const auto result = factorP * q.dot( fDirection );
        return result;
      };

      RHSisZero = true;
   }
   else
   {
      WALBERLA_ABORT( "Invalid scenario." )
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   solve( storage,
          discretization,
          solutionU,
          solutionV,
          solutionW,
          solutionP,
          initialU,
          initialV,
          initialW,
          initialP,
          rhsU,
          rhsV,
          rhsW,
          minLevel,
          maxLevel,
          multigridSettings,
          smootherSettings,
          coarseGridSettings,
          projectPressure,
          projectPressureAfterRestriction,
          vtk,
          benchmarkName,
          db,
          timingFile,
          RHSisZero );
}

} // namespace scaling_workshop
} // namespace hyteg

int main( int argc, char** argv )
{
   hyteg::scaling_workshop::benchmark( argc, argv );
}
