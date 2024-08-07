/*
 * Copyright (c) 2024 Michael Zikeli.
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

#define WRITE_BENCHMARK
#define USE_WEAK_NORM
#undef WRITE_VTK
#undef PERFORMANCE_RUN
#undef DEACTIVATE_FLOATING_CHECKS
#undef CHECK_FOR_UNDERFLOW

// +++ Operators +++
#include <ctime>

#include "Setups.h"
#include "operators-used/P1ElementwiseDiffusion_cubes_const_float16.hpp"
#include "operators-used/P1ElementwiseDiffusion_cubes_const_float32.hpp"
#include "operators-used/P1ElementwiseDiffusion_cubes_const_float64.hpp"
#include "solvePDEOpGen.hpp"

int main( int argc, char* argv[] )
{
   // +++ Setup Environment +++
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );
   walberla::MPIManager::instance()->useWorldComm();

#ifndef PERFORMANCE_RUN
   std::ofstream of( "ChebyCoeffs.dat", std::ios::trunc );
   { // Print a timestamp
      time_t     rawtime;
      struct tm* timeinfo;
      char       buffer[80];

      time( &rawtime );
      timeinfo = localtime( &rawtime );

      strftime( buffer, sizeof( buffer ), "%d-%m-%Y %H:%M:%S", timeinfo );
      std::string str( buffer );

      of << str << std::endl;
   }
#endif // n-def PERFORMANCE_RUN

   //   constexpr uint_t Quad = 0; // Order of quadrature for Diffusion (P1->0, P2->2)
   // The only Orders for quadrature that are implemented are Order P1->5 and P2->7
   // Since I'm no finite element expert and am not so sure what how to add the other orders, I'll just use these two since it's just the initialization.
   constexpr uint_t Quad = 5;

   const std::array< std::tuple< uint_t, double >, 3 > levelDependentParameters = {
       std::make_tuple( uint_t( 6 ), static_cast< double >( 1.5 * 6.206666e-06 ) ),
       std::make_tuple( uint_t( 7 ), static_cast< double >( 1.5 * 1.540316e-06 ) ),
       std::make_tuple( uint_t( 8 ), static_cast< double >( 1.5 * 3.840201e-07 ) ) };
   using SetupType =
       hyteg::mixedPrecisionMT::Setup< hyteg::mixedPrecisionMT::PoissonSetup, hyteg::mixedPrecisionMT::SolverRegister::GMG >;

   for ( const auto& [investigationLevel, errorThreshold] : levelDependentParameters )
   {
      const SetupType config = {
          .minLevel                  = uint_t( investigationLevel ),
          .maxLevel                  = uint_t( investigationLevel ),
          .errorCalcLevels           = uint_t( 1 ),
          .maxOuterSolverSteps       = uint_t( 100 ),
          .maxUnnecessaryRefinements = uint_t( 0 ), // number of Refinements after the disc.-error does not decrease anymore
          .errorThreshold            = walberla::real_t( errorThreshold ),
          .stoppingCriteria =
              hyteg::mixedPrecisionMT::setupMasks::ErrorThreshold | hyteg::mixedPrecisionMT::setupMasks::Iterations,
          .residualNormType = hyteg::mixedPrecisionMT::setupMasks::Relative | hyteg::mixedPrecisionMT::setupMasks::Vector,
          .errorNormType = hyteg::mixedPrecisionMT::setupMasks::Absolute | hyteg::mixedPrecisionMT::setupMasks::NormalizedVector,
          // +++ Solver setups +++
          .numerOfInnerSolveCyclesIR = uint_t( 1 ),
          .preSmoothingStepsGMG      = uint_t( 1 ),
          .postSmoothingStepsGMG     = uint_t( 1 ),
          .coarseLevelGMG            = uint_t( 1 ),
          .polynomialOrderCheby      = uint_t( 1 ),
          .cgIterations              = uint_t( 5 ),
          .problemSetup =
              {
                  .dim3 = true,
              },
      };
      config.setupCorrectnessCheck();

      // +++ Use Float32 as inner Solve precision +++
      hyteg::mixedPrecisionMT::solveLSEWithOperatorFold<
          SetupType,
          Quad,
          hyteg::mixedPrecisionMT::SolverRegister::IR,
          /* Smoother Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float32,
          // ++++++++++++++++++++++++ //
          /* Residual Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float64 >( config );

      // +++ Use Float32 precision GMG Solver +++
      hyteg::mixedPrecisionMT::solveLSEWithOperatorFold<
          SetupType,
          Quad,
          hyteg::mixedPrecisionMT::SolverRegister::GMG,
          /* Smoother Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float32,
          // ++++++++++++++++++++++++ //
          /* Residual Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float32 >( config );

      // +++ Use Float64 precision GMG Solver +++
      hyteg::mixedPrecisionMT::solveLSEWithOperatorFold<
          SetupType,
          Quad,
          hyteg::mixedPrecisionMT::SolverRegister::GMG,
          /* Smoother Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float64,
          // ++++++++++++++++++++++++ //
          /* Residual Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float64 >( config );
   }

   return EXIT_SUCCESS;
}
