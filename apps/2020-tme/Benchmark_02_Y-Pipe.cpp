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

using walberla::math::pi;

std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > exactW = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& ) {
   return 0;
};

std::function< real_t( const hyteg::Point3D& ) > rhsU = []( const hyteg::Point3D& ) {
   return 0;
};
std::function< real_t( const hyteg::Point3D& ) > rhsV = []( const hyteg::Point3D& ) {
   return 0;
};
std::function< real_t( const hyteg::Point3D& ) > rhsW = []( const hyteg::Point3D& ) {
   return 0;
};

void benchmark( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petScManager( &argc, &argv );

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./Benchmark_02_Y-Pipe.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );
   const walberla::Config::BlockHandle discrErrorConf = mainConf.getBlock( "DiscretizationErrorSolver" );

   const std::string discretizationString = mainConf.getParameter< std::string >( "discretization" );

   const uint_t minLevel        = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel        = mainConf.getParameter< uint_t >( "maxLevel" );

   const bool        vtk    = mainConf.getParameter< bool >( "vtk" );
   const std::string dbFile = mainConf.getParameter< std::string >( "dbFile" );

   const uint_t normCalculationLevelIncrement = mainConf.getParameter< uint_t >( "normCalculationLevelIncrement" );

   const bool solveWithCoarseGridSolverOnEachFMGLevel = mainConf.getParameter< bool >( "solveWithCoarseGridSolverOnEachFMGLevel" );

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

   const bool calculateDiscretizationError = mainConf.getParameter< bool >( "calculateDiscretizationError" );

   MultigridSettings multigridSettingsDiscrError;
   multigridSettingsDiscrError.preSmooth                 = discrErrorConf.getParameter< uint_t >( "preSmooth" );
   multigridSettingsDiscrError.postSmooth                = discrErrorConf.getParameter< uint_t >( "postSmooth" );
   multigridSettingsDiscrError.incSmooth                 = discrErrorConf.getParameter< uint_t >( "incSmooth" );
   multigridSettingsDiscrError.fmgInnerIterations        = discrErrorConf.getParameter< uint_t >( "fmgInnerIterations" );
   multigridSettingsDiscrError.numCycles                 = discrErrorConf.getParameter< uint_t >( "numCycles" );
   multigridSettingsDiscrError.absoluteResidualTolerance = discrErrorConf.getParameter< real_t >( "absoluteResidualTolerance" );

   SmootherSettings smootherSettingsDiscrError;
   smootherSettingsDiscrError.estimateOmega             = discrErrorConf.getParameter< bool >( "estimateOmega" );
   smootherSettingsDiscrError.omega                     = discrErrorConf.getParameter< real_t >( "omega" );
   smootherSettingsDiscrError.omegaEstimationLevel      = discrErrorConf.getParameter< uint_t >( "omegaEstimationLevel" );
   smootherSettingsDiscrError.omegaEstimationIterations = discrErrorConf.getParameter< uint_t >( "omegaEstimationIterations" );
   smootherSettingsDiscrError.numGSVelocity             = discrErrorConf.getParameter< uint_t >( "numGSVelocity" );
   smootherSettingsDiscrError.symmGSVelocity            = discrErrorConf.getParameter< bool >( "symmGSVelocity" );

   CoarseGridSettings coarseGridSettingsDiscrError;
   coarseGridSettingsDiscrError.absoluteResidualTolerance = discrErrorConf.getParameter< real_t >( "coarseGridAbsoluteResidualTolerance" );
   coarseGridSettingsDiscrError.maxIterations             = discrErrorConf.getParameter< uint_t >( "maxIterations" );
   coarseGridSettingsDiscrError.solverType                = discrErrorConf.getParameter< uint_t >( "coarseGridSolverType" );

   const std::string dbFileDiscrError = discrErrorConf.getParameter< std::string >( "dbFile" );
   

   Discretization discretization = Discretization::P2_P1;
   if ( discretizationString == "p1p1" )
   {
      discretization = Discretization::P1_P1;
   }

   int tWidth = 6;
   int tBase = 4;
   const real_t alpha = (0.5/real_c(tWidth-1)) * pi;
   real_t eps = 1e-3;
   real_t tiltZSlope = 0.2;


   std::set< std::array< int, 3 > > cubeCoords;
   for ( int i = 0; i < tWidth; i++ )
   {
      cubeCoords.insert( {-i, 0, 0} );
      cubeCoords.insert( {i, 0, 0} );
   }
   for ( int i = 1; i < tBase; i++ )
   {
      cubeCoords.insert( {0, -i, 0} );
   }

   auto meshInfo = MeshInfo::meshCubedDomain( cubeCoords, 1 );

   auto outflowBoundary = [&]( const Point3D & p )
   {
     if ( p[0] > real_c( tWidth ) - eps && p[1] > -eps && std::abs( p[2] - 0.5 ) < eps )
     {
        return true;
     }
     else if ( p[0] < - real_c( tWidth ) + 1 + eps && p[1] > -eps && std::abs( p[2] - 0.5 ) < eps )
     {
        return true;
     }

     return false;
   };

   auto tiltZCoordinateMap = [&]( const hyteg::Point3D& p ) -> Point3D {

     if ( p[0] > 1 + eps && p[1] > -eps )
     {
        auto pNew = p;
        pNew[2] = p[2] + tiltZSlope * (p[0]-1);
        return pNew;
     }

     if ( p[0] < -eps && p[1] > -eps )
     {
        auto pNew = p;
        pNew[2] = p[2] + tiltZSlope * p[0];
        return pNew;
     }

     return p;
   };

   meshInfo.applyCoordinateMap( tiltZCoordinateMap );

   real_t tDist = real_c( tWidth - 1 );

   auto pipeCoordinateMap = [&]( const hyteg::Point3D& p ) -> Point3D {

      if ( p[0] > tDist - eps && p[1] > -eps )
      {
         // right anchor
         Point3D anchorRight( {tDist, 0, 0} );
         auto     pShift = p - anchorRight;
         Matrix3r m( Matrix3r::Zero() );
         m(0, 0) = std::cos(alpha);
         m(0, 1) = -std::sin(alpha);
         m(1, 0) = std::sin(alpha);
         m(1, 1) = std::cos(alpha);
         m(2, 2) = 1;
         auto pShiftRotated = m.mul(pShift);
         if ( pShiftRotated[1] < pShift[1] )
         {
            pShiftRotated[1] = pShift[1];
         }
         auto pShiftRotatedBackshift = pShiftRotated + anchorRight;
         return pShiftRotatedBackshift;
      }

     if ( p[0] < eps - tDist + 1 && p[1] > -eps )
     {
        // left anchor
        Point3D anchorLeft( {-tDist + 1, 0, 0} );
        auto     pShift = p - anchorLeft;
        Matrix3r m( Matrix3r::Zero() );
        m(0, 0) = std::cos(-alpha);
        m(0, 1) = -std::sin(-alpha);
        m(1, 0) = std::sin(-alpha);
        m(1, 1) = std::cos(-alpha);
        m(2, 2) = 1;
        auto pShiftRotated = m.mul(pShift);
        if ( pShiftRotated[1] < pShift[1] )
        {
           pShiftRotated[1] = pShift[1];
        }
        auto pShiftRotatedBackshift = pShiftRotated + anchorLeft;
        return pShiftRotatedBackshift;
     }

      return p;
   };

   for ( int i = 0; i < tWidth-1; i++ )
   {
      tDist = real_c(tWidth - i);
      meshInfo.applyCoordinateMap( pipeCoordinateMap );
   }

   std::function< real_t( const hyteg::Point3D& ) > exactV = [&]( const hyteg::Point3D& p ) {
     if ( std::abs( p[1] - real_c( -tBase + 1 ) ) < eps && p[0] > -eps && p[0] < 1 + eps )
     {
        const Point3D center( {0.5, real_c( -tBase + 1 ), 0.5} );
        const auto    radius  = 0.5;
        const auto    shifted = ( p - center ) / radius;

        return ( 1 - std::sin( 0.5 * pi * shifted[0] * shifted[0] ) ) *
               ( 1 - std::sin( 0.5 * pi * shifted[2] * shifted[2] ) );
     }
     else
     {
        return 0.0;
     }
   };

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // new code ...
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( 2, outflowBoundary, false );

   // ... replaces old
   // meshInfo.setAllMeshBoundaryFlags( 1 );
   // meshInfo.setMeshBoundaryFlagsByVertexLocation( 2, outflowBoundary, false );
   // setupStorage->setMeshBoundaryFlagsInner( 0, true );

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
          multigridSettingsDiscrError,
          smootherSettingsDiscrError,
          coarseGridSettingsDiscrError,
          false,
          false,
          calculateDiscretizationError,
          normCalculationLevelIncrement,
          solveWithCoarseGridSolverOnEachFMGLevel,
          vtk,
          "Benchmark_02_Y-Pipe",
          false,
          dbFile,
          dbFileDiscrError );
}

} // namespace tme_benchmarks
} // namespace hyteg

int main( int argc, char** argv )
{
   hyteg::tme_benchmarks::benchmark( argc, argv );
}
