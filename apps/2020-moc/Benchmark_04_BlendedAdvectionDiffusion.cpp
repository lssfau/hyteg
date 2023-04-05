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
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "Helpers.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

namespace hyteg {
namespace moc_benchmarks {

const real_t INITIAL_DIFFUSIVITY_TIME_PRODUCT = 1e-03 * 2 * pi;

class TempSolution : public Solution
{
 public:
   TempSolution( real_t diffusivity, Point3D p0, real_t t0 )
   : Solution( t0 )
   , diffusivity_( diffusivity )
   , p0_( p0 )
   {}

   real_t operator()( const Point3D& x ) const override
   {
      auto x_hat = p0_[0] * std::cos( currentTime_ ) - p0_[1] * std::sin( currentTime_ );
      auto y_hat = -p0_[0] * std::sin( currentTime_ ) + p0_[1] * std::cos( currentTime_ );
      //      auto x_hat    = p0_[0] + currentTime_;
      //      auto y_hat    = 0;
      auto exponent = -std::pow( r( x, x_hat, y_hat ), 2 ) / ( 4.0 * diffusivity_ * currentTime_ );
      return ( 1.0 / ( 4.0 * pi * diffusivity_ * currentTime_ ) ) * std::exp( exponent );
   }

 private:
   real_t r( const Point3D& p, const real_t& x_hat, const real_t& y_hat ) const
   {
      return std::sqrt( std::pow( p[0] - x_hat, 2 ) + std::pow( p[1] - y_hat, 2 ) + std::pow( p[2], 2 ) );
   }

 private:
   real_t  diffusivity_;
   Point3D p0_;
};

class VelocitySolutionX : public Solution
{
   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& x ) const override { return -x[1]; }
};

class VelocitySolutionY : public Solution
{
   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& x ) const override { return x[0]; }
};

class VelocitySolutionZ : public Solution
{
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
      auto defaultFile = "./Benchmark_04_BlendedAdvectionDiffusion.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t      numTimeSteps      = mainConf.getParameter< uint_t >( "numTimeSteps" );
   const bool        threeDim          = mainConf.getParameter< bool >( "threeDim" );
   const uint_t      level             = mainConf.getParameter< uint_t >( "level" );
   const real_t      diffusivity       = mainConf.getParameter< real_t >( "diffusivity" );
   const bool        rotationOnly      = mainConf.getParameter< bool >( "rotationOnly" );
   const bool        resetParticles    = mainConf.getParameter< bool >( "resetParticles" );
   const bool        adjustedAdvection = mainConf.getParameter< bool >( "adjustedAdvection" );
   const bool        crankNicolson     = mainConf.getParameter< bool >( "crankNicolson" );
   const bool        strangSplitting   = mainConf.getParameter< bool >( "strangSplitting" );
   const uint_t      printInterval     = mainConf.getParameter< uint_t >( "printInterval" );
   const bool        vtk               = mainConf.getParameter< bool >( "vtk" );
   const uint_t      vtkInterval       = mainConf.getParameter< uint_t >( "vtkInterval" );
   const bool        verbose           = mainConf.getParameter< bool >( "verbose" );
   const std::string dbFile            = mainConf.getParameter< std::string >( "dbFile" );

   LoadBalancingOptions lbOptions;
   lbOptions.type = 0;

   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
   if ( threeDim )
   {
      meshInfo = MeshInfo::meshSphericalShell( 3, 3, 0.5, 1.5 );
   }
   else
   {
      meshInfo = MeshInfo::meshAnnulus( 0.5, 1.5, MeshInfo::CROSS, 6, 2 );
   }

   const Point3D p0( 0, 1, 0 );

   const real_t tStart = INITIAL_DIFFUSIVITY_TIME_PRODUCT / diffusivity;
   const real_t tEnd   = tStart + 2.0 * pi;

   const real_t dt = ( tEnd - tStart ) / real_c( numTimeSteps );

   TempSolution      cSolution( diffusivity, p0, tStart );
   VelocitySolutionX uSolution;
   VelocitySolutionY vSolution;
   VelocitySolutionZ wSolution;

   solve< P2Function< real_t >,
          P2ElementwiseBlendingLaplaceOperator,
          P2ElementwiseBlendingMassOperator,
          P2ElementwiseUnsteadyDiffusionOperator >( meshInfo,
                                                    true,
                                                    cSolution,
                                                    uSolution,
                                                    vSolution,
                                                    wSolution,
                                                    dt,
                                                    diffusivity,
                                                    level,
                                                    crankNicolson ? DiffusionTimeIntegrator::CrankNicolson :
                                                                    DiffusionTimeIntegrator::ImplicitEuler,
                                                    !rotationOnly,
                                                    strangSplitting,
                                                    !rotationOnly || resetParticles,
                                                    1,
                                                    adjustedAdvection,
                                                    false,
                                                    false,
                                                    numTimeSteps,
                                                    lbOptions,
                                                    vtk,
                                                    true,
                                                    "Benchmark_04_BlendedAdvectionDiffusion",
                                                    printInterval,
                                                    vtkInterval,
                                                    verbose,
                                                    dbFile );
}
} // namespace moc_benchmarks
} // namespace hyteg

int main( int argc, char* argv[] )
{
   hyteg::moc_benchmarks::benchmark( argc, argv );
   return EXIT_SUCCESS;
}