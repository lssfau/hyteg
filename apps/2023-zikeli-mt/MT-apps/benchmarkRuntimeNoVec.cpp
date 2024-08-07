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

#define PERFORMANCE_RUN
#define DEACTIVATE_FLOATING_CHECKS
#define WRITE_BENCHMARK
#undef WRITE_VTK
#undef CHECK_FOR_UNDERFLOW
#undef USE_WEAK_NORM

// +++ Operators +++
#include "Setups.h"
#include "operators-used/P1ElementwiseDiffusion_cubes_const_float32.hpp"
#include "operators-used/P1ElementwiseDiffusion_cubes_const_float64.hpp"
#include "solvePDEOpGen.hpp"

int main( int argc, char* argv[] )
{
   // +++ Setup Environment +++
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );
   walberla::MPIManager::instance()->useWorldComm();

   using SetupType =
       hyteg::mixedPrecisionMT::Setup< hyteg::mixedPrecisionMT::PoissonSetup, hyteg::mixedPrecisionMT::SolverRegister::GMG >;
   // level |  errorThreshold
   //      const std::array<std::tuple<uint_t, double>, 4> levelDependentParameters = { std::make_tuple( uint_t(6), static_cast<double>(1.5*6.206666e-06) ),
   //                                                                                   std::make_tuple( uint_t(7), static_cast<double>(1.5*1.540316e-06) ),
   //                                                                                   std::make_tuple( uint_t(8), static_cast<double>(1.5*3.840201e-07) ),
   //                                                                                   std::make_tuple( uint_t(9), static_cast<double>(1.5*9.594706e-08) ) };
   const std::array< std::tuple< uint_t, double >, 3 > levelDependentParameters = {
       std::make_tuple( uint_t( 6 ), static_cast< double >( 1.5 * 6.206666e-06 ) ),
       std::make_tuple( uint_t( 7 ), static_cast< double >( 1.5 * 1.540316e-06 ) ),
       std::make_tuple( uint_t( 8 ), static_cast< double >( 1.5 * 3.840201e-07 ) ) };
   //   const std::array<std::tuple<uint_t, double>, 1> levelDependentParameters = { std::make_tuple( uint_t(9), static_cast<double>(1.5*9.594706e-08) ) };
   const std::array< uint_t, 2 > numberVCycles = { 1, 4 };
   for ( const auto& [investigationLevel, errorThreshold] : levelDependentParameters )
   {
      for ( uint_t numberVCycle : numberVCycles )
      {
         const SetupType config = {
             .minLevel                  = uint_t( investigationLevel ),
             .maxLevel                  = uint_t( investigationLevel ),
             .maxOuterSolverSteps       = uint_t( 100 ),
             .maxUnnecessaryRefinements = uint_t( 0 ), // number of Refinements after the disc.-error does not decrease anymore
             .errorThreshold            = walberla::real_t( errorThreshold ),
             .stoppingCriteria =
                 hyteg::mixedPrecisionMT::setupMasks::ErrorThreshold | hyteg::mixedPrecisionMT::setupMasks::Iterations,
             .residualNormType = hyteg::mixedPrecisionMT::setupMasks::Relative | hyteg::mixedPrecisionMT::setupMasks::Vector,
             .errorNormType =
                 hyteg::mixedPrecisionMT::setupMasks::Absolute | hyteg::mixedPrecisionMT::setupMasks::NormalizedVector,
             // +++ Solver setups +++
             .numerOfInnerSolveCyclesIR = uint_t( numberVCycle ),
             .preSmoothingStepsGMG      = uint_t( 1 ),
             .postSmoothingStepsGMG     = uint_t( 1 ),
             .coarseLevelGMG            = uint_t( 1 ),
             .polynomialOrderCheby      = uint_t( 1 ),
             .cgIterations              = uint_t( 5 ),
             // +++ Output Names +++
             .benchFile = "Runtime",
             .problemSetup =
                 {
                     .dim3 = true,
                 },
         };
         config.setupCorrectnessCheck();

         //   constexpr uint_t Quad = 0; // Order of quadrature for Diffusion (P1->0, P2->2)
         // The only Orders for quadrature that are implemented are Order P1->5 and P2->7
         // Since I'm no finite element expert and am not so sure what how to add the other orders, I'll just use these two since it's just the initialization.
         constexpr uint_t Quad = 5;

         if ( numberVCycle == 1 )
         {
            if ( investigationLevel <= 7 )
            {
               // +++ Use Float32 precision GMG Solver +++
               hyteg::mixedPrecisionMT::solveLSEWithOperatorFold<
                   SetupType,
                   Quad,
                   hyteg::mixedPrecisionMT::SolverRegister::GMG,
                   /* Smoother Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float32,
                   // ++++++++++++++++++++++++ //
                   /* Residual Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float32 >(
                   config );
            } // if ( investigationLevel <= 7 )

            // +++ Use Float64 precision GMG Solver +++
            hyteg::mixedPrecisionMT::solveLSEWithOperatorFold<
                SetupType,
                Quad,
                hyteg::mixedPrecisionMT::SolverRegister::GMG,
                /* Smoother Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float64,
                // ++++++++++++++++++++++++ //
                /* Residual Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float64 >( config );
         } // if ( innerVCycle == 1 )

         // +++ Use Float32 as inner Solve precision +++
         hyteg::mixedPrecisionMT::solveLSEWithOperatorFold<
             SetupType,
             Quad,
             hyteg::mixedPrecisionMT::SolverRegister::IR,
             /* Smoother Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float32,
             // ++++++++++++++++++++++++ //
             /* Residual Operator for IR */ hyteg::operatorgeneration::P1ElementwiseDiffusion_cubes_const_float64 >( config );
      } // for ( v-cycles )
   }    // for ( investigationLevel )

   return EXIT_SUCCESS;
}