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
#include "constantStencilOperator/P2ConstantOperator.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "constantStencilOperator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {
namespace moc_benchmarks {

auto r = []( const hyteg::Point3D& x, const hyteg::Point3D& x0, const real_t& r0 ) -> real_t {
   return ( 1 / r0 ) * std::sqrt( std::pow( x[0] - x0[0], 2 ) + std::pow( x[1] - x0[1], 2 ) + std::pow( x[2] - x0[2], 2 ) );
};

auto conicalBody = []( const hyteg::Point3D& x ) -> real_t {
   const Point3D x0( 0.5, 0.25, 0.0 );
   const real_t  r0 = 0.15;
   if ( r( x, x0, r0 ) <= 1. )
      return 1 - r( x, x0, r0 );
   else
      return 0.0;
};

auto gaussianCone = []( const hyteg::Point3D& x ) -> real_t {
   const Point3D x0( 0.25, 0.5, 0.0 );
   const real_t  r0 = 0.15;
   if ( r( x, x0, r0 ) <= 1. )
      return ( 1 + std::cos( walberla::math::pi * r( x, x0, r0 ) ) ) * 0.25;
   else
      return 0.0;
};

auto slottedCylinder = []( const hyteg::Point3D& x ) -> real_t {
   const Point3D x0( 0.5, 0.75, 0.0 );
   const real_t  r0 = 0.15;
   if ( ( r( x, x0, r0 ) <= 1. ) && ( std::abs( x[0] - x0[0] ) >= 0.025 || x[1] >= 0.85 ) )
      return 1;
   else
      return 0.0;
};

class TempSolution : public Solution
{
 public:
   TempSolution( bool enableGaussianCone, bool enableLinearCone, bool enableCylinder )
   : enableGaussianCone_( enableGaussianCone )
   , enableLinearCone_( enableLinearCone )
   , enableCylinder_( enableCylinder )
   {}

   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& x ) const override
   {
      real_t val = 0;
      if ( enableGaussianCone_ )
         val += gaussianCone( x );
      if ( enableLinearCone_ )
         val += conicalBody( x );
      if ( enableCylinder_ )
         val += slottedCylinder( x );
      return val;
   }

 private:
   bool enableGaussianCone_;
   bool enableLinearCone_;
   bool enableCylinder_;
};

class VelocitySolutionX : public Solution
{
   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& x ) const override { return 0.5 - x[1]; }
};

class VelocitySolutionY : public Solution
{
   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& x ) const override { return x[0] - 0.5; }
};

class VelocitySolutionZ : public Solution
{
 public:
   explicit VelocitySolutionZ( bool threeDim )
   : threeDim_( threeDim )
   {}

   /// Evaluates the solution at a specific point.
   real_t operator()( const Point3D& x ) const override
   {
      if ( !threeDim_ )
      {
         return 0;
      }

      const auto xx = x[0] - 0.5;
      return 0.1 * xx;
   }

 private:
   bool threeDim_;
};

void benchmark( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./Benchmark_01_CircularAdvection.prm";
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
   const bool        threeDim               = mainConf.getParameter< bool >( "threeDim" );
   const bool        resetParticles         = mainConf.getParameter< bool >( "resetParticles" );
   const uint_t      resetParticlesInterval = mainConf.getParameter< uint_t >( "resetParticlesInterval" );
   const bool        adjustedAdvection      = mainConf.getParameter< bool >( "adjustedAdvection" );
   const bool        enableCylinder         = mainConf.getParameter< bool >( "enableCylinder" );
   const bool        enableLinearCone       = mainConf.getParameter< bool >( "enableLinearCone" );
   const bool        enableGaussianCone     = mainConf.getParameter< bool >( "enableGaussianCone" );
   const uint_t      printInterval          = mainConf.getParameter< uint_t >( "printInterval" );
   const bool        vtk                    = mainConf.getParameter< bool >( "vtk" );
   const uint_t      vtkInterval            = mainConf.getParameter< uint_t >( "vtkInterval" );
   const std::string dbFile                 = mainConf.getParameter< std::string >( "dbFile" );
   const bool        globalMaxLimiter       = mainConf.getParameter< bool >( "globalMaxLimiter" );
   const std::string spaceDiscretization    = mainConf.getParameter< std::string >( "spaceDiscretization" );

   LoadBalancingOptions lbOptions;
   lbOptions.type = 0;

   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
   if ( threeDim )
   {
      meshInfo = MeshInfo::meshCuboid( Point3D( 0, 0, -0.5 ), Point3D( 1, 1, 0.5 ), 1, 1, 1 );
   }
   else
   {
      meshInfo = MeshInfo::meshRectangle( Point2D( 0, 0 ), Point2D( 1, 1 ), MeshInfo::CRISS, 1, 1 );
   }

   const real_t tEnd = 2 * walberla::math::pi;
   const real_t dt   = tEnd / real_c( numTimeSteps );

   TempSolution      cSolution( enableGaussianCone, enableLinearCone, enableCylinder );
   VelocitySolutionX uSolution;
   VelocitySolutionY vSolution;
   VelocitySolutionZ wSolution( threeDim );

   if ( spaceDiscretization == "P1" )
   {
      solve< P1Function< real_t >, P1ConstantLaplaceOperator, P1ConstantMassOperator, P1ConstantUnsteadyDiffusionOperator >(
          meshInfo,
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
          globalMaxLimiter,
          true,
          numTimeSteps,
          lbOptions,
          vtk,
          true,
          "Benchmark_01_CircularAdvection",
          printInterval,
          vtkInterval,
          false,
          dbFile );
   }
   else
   {
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
                                                       globalMaxLimiter,
                                                       true,
                                                       numTimeSteps,
                                                       lbOptions,
                                                       vtk,
                                                       true,
                                                       "Benchmark_01_CircularAdvection",
                                                       printInterval,
                                                       vtkInterval,
                                                       false,
                                                       dbFile );
   }
}
} // namespace moc_benchmarks
} // namespace hyteg

int main( int argc, char* argv[] )
{
   hyteg::moc_benchmarks::benchmark( argc, argv );
   return EXIT_SUCCESS;
}