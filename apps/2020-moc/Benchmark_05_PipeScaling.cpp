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

#include <core/Environment.h>
#include <core/math/Constants.h>

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "Helpers.hpp"
#include "constant_stencil_operator/P2ConstantOperator.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

namespace hyteg {
namespace moc_benchmarks {

auto r = []( const hyteg::Point3D& x, const hyteg::Point3D& x0, const real_t& r0 ) -> real_t {
   return ( 1 / r0 ) * std::sqrt( std::pow( x[0] - x0[0], 2 ) + std::pow( x[1] - x0[1], 2 ) + std::pow( x[2] - x0[2], 2 ) );
};

auto gaussianCone = []( const hyteg::Point3D& x ) -> real_t {
   const Point3D x0( 0.5, 0.5, 0.5 );
   const real_t  r0 = 0.15;
   if ( r( x, x0, r0 ) <= 1. )
      return ( 1 + std::cos( walberla::math::pi * r( x, x0, r0 ) ) ) * 0.25;
   else
      return 0.0;
};

class TempSolution : public Solution
{
 public:
   TempSolution() {}

   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& x ) const override
   {
      auto xTranslated = x;
      while ( xTranslated[0] > 1.0 )
      {
         xTranslated[0] -= 1.0;
      }
      return gaussianCone( xTranslated );
   }
};

class VelocitySolutionX : public Solution
{
 public:
   explicit VelocitySolutionX()
   : Solution()
   {}

   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& ) const override { return 1.0; }
};

class VelocitySolutionY : public Solution
{
 public:
   explicit VelocitySolutionY()
   : Solution()
   {}

   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& ) const override { return 0; }
};

class VelocitySolutionZ : public Solution
{
 public:
   explicit VelocitySolutionZ()
   : Solution()
   {}

   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& ) const override { return 0; }
};

void benchmark( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./Benchmark_05_PipeScaling.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t      numTimeSteps           = mainConf.getParameter< uint_t >( "numTimeSteps" );
   const uint_t      level                  = mainConf.getParameter< uint_t >( "level" );
   const bool        resetParticles         = mainConf.getParameter< bool >( "resetParticles" );
   const uint_t      resetParticlesInterval = mainConf.getParameter< uint_t >( "resetParticlesInterval" );
   const bool        adjustedAdvection      = mainConf.getParameter< bool >( "adjustedAdvection" );
   const uint_t      printInterval          = mainConf.getParameter< uint_t >( "printInterval" );
   const bool        vtk                    = mainConf.getParameter< bool >( "vtk" );
   const uint_t      vtkInterval            = mainConf.getParameter< uint_t >( "vtkInterval" );
   const std::string dbFile                 = mainConf.getParameter< std::string >( "dbFile" );
   const uint_t      diameterCubes          = mainConf.getParameter< uint_t >( "diameterCubes" );
   const uint_t      lengthCubes            = mainConf.getParameter< uint_t >( "lengthCubes" );
   const real_t      dt                     = mainConf.getParameter< real_t >( "timeStepSize" );

   LoadBalancingOptions lbOptions;
   lbOptions.type = mainConf.getParameter< uint_t >( "lbType" );

   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   WALBERLA_CHECK_EQUAL( lengthCubes % diameterCubes, 0 );

   meshInfo = MeshInfo::meshCuboid( Point3D( 0, 0, 0 ),
                                    Point3D( real_c( lengthCubes / diameterCubes ), 1, 1 ),
                                    lengthCubes,
                                    diameterCubes,
                                    diameterCubes );

   TempSolution      cSolution;
   VelocitySolutionX uSolution;
   VelocitySolutionY vSolution;
   VelocitySolutionZ wSolution;

   WALBERLA_LOG_INFO_ON_ROOT( " - length in cubes (x):      " << lengthCubes );
   WALBERLA_LOG_INFO_ON_ROOT( " - diameter in cubes (y, z): " << diameterCubes );
   WALBERLA_LOG_INFO_ON_ROOT( " - LB type:                  " << lbOptions.type );

   solve< P2Function< real_t >,
          P2ElementwiseBlendingLaplaceOperator,
          P2ElementwiseBlendingMassOperator,
          P2ElementwiseUnsteadyDiffusionOperator >( meshInfo,
                                                    false,
                                                    cSolution,
                                                    uSolution,
                                                    vSolution,
                                                    wSolution,
                                                    dt,
                                                    1.0,
                                                    level,
                                                    DiffusionTimeIntegrator::ImplicitEuler,
                                                    false,
                                                    false,
                                                    resetParticles,
                                                    resetParticlesInterval,
                                                    adjustedAdvection,
                                                    false,
                                                    false,
                                                    numTimeSteps,
                                                    lbOptions,
                                                    vtk,
                                                    true,
                                                    "Benchmark_05_PipeScaling",
                                                    printInterval,
                                                    vtkInterval,
                                                    true,
                                                    dbFile );
}
} // namespace moc_benchmarks
} // namespace hyteg

int main( int argc, char* argv[] )
{
   hyteg::moc_benchmarks::benchmark( argc, argv );
   return EXIT_SUCCESS;
}