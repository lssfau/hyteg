/*
* Copyright (c) 2024 Michael Zikeli.
*
* This file is part of HyTeG
* (see https://i10git.cs.fau.de/hyteg/hyteg).
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by * the Free Software Foundation, either version 3 of the License, or
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

#pragma once

#include "core/debug/CheckFunctions.h"

#include "SolverRegister.h"
#include "globalHeader.h"

namespace hyteg::mixedPrecisionMT {

using hyteg::Point3D;
using walberla::numeric_cast;

struct ProblemSetup
{
   // +++ Parameters +++
   const std::string fileDescription = "GenericProblem";
   const std::string descriptionPDE  = "GenericProblem";

   // +++ Mesh setups +++
   std::function< std::string( const bool ) > chooseMesh = []( const bool use3D ) {
      // TODO use another mesh according to the mesh generation functions.
      //      return (use3D) ? "../../../data/meshes/3D/cube_6el.msh"
      //                     : "../../../data/meshes/quad_2el.msh";
      return ( use3D ) ? "../../../data/meshes/3D/cube_24el.msh" : "../../../data/meshes/quad_8el.msh";
   };
   const bool        dim3;
   const std::string meshFile = chooseMesh( dim3 );

   // +++ Initialization Functions +++
   const std::function< real_t( const Point3D& ) > initExact = []( const Point3D& ) {
      WALBERLA_ABORT( "There is no exact solution for a generic setup." );
      return real_c( 0.0 );
   };

   const std::function< real_t( const Point3D& ) > initRHS = []( const Point3D& ) {
      WALBERLA_ABORT( "There is no right hand side solution for a generic setup." );
      return real_c( 0.0 );
   };
};

struct LaplaceSetup : ProblemSetup
{
   // +++ Parameters +++
   const std::string fileDescription = "SolveLaplace";
   const std::string descriptionPDE  = "Laplace-UnitCube-RHS=0-u=sin(x)sinh(y)z";

   // +++ Mesh setups +++
   const bool        dim3     = ProblemSetup::dim3;
   const std::string meshFile = ProblemSetup::chooseMesh( dim3 );

   // +++ Initialization Functions +++
   const std::function< real_t( const Point3D& ) > initExact = [this]( const Point3D& p ) {
      if ( this->dim3 )
      {
         return real_c( sin( p[0] ) * sinh( p[1] ) * p[2] );
      }
      else
      {
         return real_c( sin( p[0] ) * sinh( p[1] ) );
      }
   };

   const std::function< real_t( const Point3D& ) > initRHS = []( const Point3D& ) {
      return real_c( 0.0 ); // ( laplace( sin(x)sinh(y)z ) = 0 // Wolfram Alpha )
   };

}; // struct LaplaceSetup

struct PoissonSetup : ProblemSetup
{
   // +++ Parameters +++
   const std::string fileDescription = "SolvePoisson";
   const std::string descriptionPDE  = "Poisson-UnitCube-RHS=sin(2x)sinh(y)cosh(z)-u=(-0.5sin(2x)sinh(y)cosh(z))";

   // +++ Mesh setups +++
   const bool        dim3     = ProblemSetup::dim3;
   const std::string meshFile = ProblemSetup::chooseMesh( dim3 );

   // +++ Initialization Functions +++
   const std::function< real_t( const Point3D& ) > initExact = [this]( const Point3D& p ) {
      if ( this->dim3 )
      {
         return real_c( 0.5 * sin( 2 * p[0] ) * sinh( p[1] ) * cosh( p[2] ) );
      }
      else
      {
         return real_c( 0.5 * sin( 2 * p[0] ) * sinh( p[1] ) );
      }
   };

   const std::function< real_t( const Point3D& ) > initRHS = [this]( const Point3D& p ) {
      if ( this->dim3 )
      {
         return real_c( sin( 2 * p[0] ) * sinh( p[1] ) * cosh( p[2] ) );
      }
      else
      {
         return real_c( sin( 2 * p[0] ) * sinh( p[1] ) * 1.5 );
      }
   };

}; // struct LaplaceSetup

template < class Problem_t, SolverRegister IRSmoother_T >
struct Setup
{
   static constexpr SolverRegister IRSmoother = IRSmoother_T;

   // +++ Parameters +++
   // +++ Runtime setups +++
   const uint_t minLevel                  = uint_t( 1 );
   const uint_t maxLevel                  = uint_t( 4 ); // in 3D 8 highest possible on my notebook
   const uint_t errorCalcLevels           = uint_t( 1 );
   const uint_t maxOuterSolverSteps       = uint_t( 100 );
   const uint_t maxUnnecessaryRefinements = uint_t( 2 ); // number of Refinements after the disc.-error does not decrease anymore
   const real_t minResidualRate           = real_t( 0.91 );
   const real_t minErrorRate              = real_t( 0.91 );
   const real_t residualThreshold         = real_t( std::numeric_limits< real_t >::min() );
   const real_t errorThreshold            = real_t( std::numeric_limits< real_t >::epsilon() );
   const uint_t stoppingCriteria          = setupMasks::All;
   const uint_t residualNormType          = setupMasks::Relative | setupMasks::Vector;
   const uint_t errorNormType             = setupMasks::Absolute | setupMasks::NormalizedVector;

   // +++ Solver setups +++
   const uint_t numerOfInnerSolveCyclesIR = uint_t( 2 );
   const uint_t preSmoothingStepsGMG      = uint_t( 2 );
   const uint_t postSmoothingStepsGMG     = uint_t( 2 );
   const uint_t coarseLevelGMG            = uint_t( 1 );
   const uint_t polynomialOrderCheby      = uint_t( 1 );
   const uint_t cgIterations =
       uint_t( 5 ); // Even for float32, with 5 residual is already one order of magnitude below discretization error.

   // +++ Output setups +++
   const std::string benchPath  = "./benchmark";
   const std::string benchFile  = "Table";
   const std::string timingPath = "./output";
   const std::string configDescriptor =
       ( std::stringstream( "" ) << "_minLvl-" << minLevel << "_maxLvl-" << maxLevel << "_coarseLvl-" << coarseLevelGMG
                                 << "_innerSolves-" << numerOfInnerSolveCyclesIR << "_preSmooths-" << preSmoothingStepsGMG
                                 << "_postSmooths-" << postSmoothingStepsGMG << "_ChebyOrder-" << polynomialOrderCheby
                                 << "_CGIterations-" << cgIterations )
           .str();
   const bool printVTK = true;

   // +++ Setup functionality +++
   constexpr void setupCorrectnessCheck() const
   {
      using namespace setupMasks;
      WALBERLA_CHECK( ( residualNormType & Absolute ) ^ ( residualNormType & Relative ) );
      WALBERLA_CHECK( ( residualNormType & Vector ) ^ ( residualNormType & NormalizedVector ) ^ ( residualNormType & Weak ) );

      WALBERLA_CHECK( ( errorNormType & Absolute ) ^ ( errorNormType & Relative ) );
      WALBERLA_CHECK( ( errorNormType & Vector ) ^ ( errorNormType & NormalizedVector ) ^ ( errorNormType & Weak ) );
   }

   template < typename ResidualType, typename RateType, typename ErrorType >
   bool needFurtherOuterSolveIterations( const ResidualType& residualNorm,
                                         const RateType&     residualRate,
                                         const ErrorType&    errorNorm,
                                         const RateType&     errorRate,
                                         const uint_t        iterations,
                                         uint_t*             criteriaContainer = nullptr ) const
   {
      using namespace setupMasks;
      uint_t criteria = 0;

      if ( stoppingCriteria & Iterations && iterations > maxOuterSolverSteps )
      {
         criteria |= Iterations;
         WALBERLA_LOG_INFO_ON_ROOT( "\nSolver Terminates after " << iterations << " iterations, due to:\tMax Iterations ("
                                                                 << maxOuterSolverSteps << ") Exceeded .\n\n" );
      }
      if ( stoppingCriteria & ResidualThreshold && ( residualNorm < residualThreshold ) )
      {
         criteria |= ResidualThreshold;
         WALBERLA_LOG_INFO_ON_ROOT( "\nSolver Terminates after " << iterations << " iterations, due to:\tResidual small enough ( "
                                                                 << residualNorm << " < " << residualThreshold << " ).\n\n" );
      }
      if ( stoppingCriteria & ResidualRate && residualRate >= minResidualRate )
      {
         criteria |= ResidualRate;
         WALBERLA_LOG_INFO_ON_ROOT( "\nSolver Terminates after "
                                    << iterations << " iterations, due to:\tResidual converged ( Residual Rate = " << residualRate
                                    << " ).\n\n" );
      }
      if ( stoppingCriteria & ErrorThreshold && ( errorNorm < errorThreshold ) )
      {
         criteria |= ErrorThreshold;
         WALBERLA_LOG_INFO_ON_ROOT( "\nSolver Terminates after " << iterations << " iterations, due to:\tError small enough ( "
                                                                 << errorNorm << " < " << errorThreshold << " ).\n\n" );
      }
      if ( stoppingCriteria & ErrorRate && errorRate >= minErrorRate )
      {
         criteria |= ErrorRate;
         WALBERLA_LOG_INFO_ON_ROOT( "\nSolver Terminates after "
                                    << iterations << " iterations, due to:\tError converged ( Error Rate = " << errorRate
                                    << " ).\n\n" );
      }
      if ( std::isnan( residualNorm ) )
      {
         criteria |= ResIsNaN;
         WALBERLA_LOG_INFO_ON_ROOT( "\033[33m"
                                    << "\nSolver Terminates after " << iterations
                                    << " iterations, due to:\tResidual Norm is NaN.\n\n"
                                    << "\033[0m" );
      }
      if ( std::isnan( errorNorm ) )
      {
         criteria |= ErrIsNaN;
         WALBERLA_LOG_INFO_ON_ROOT( "\033[33m"
                                    << "\nSolver Terminates after " << iterations
                                    << " iterations, due to:\tError Norm is NaN.\n\n"
                                    << "\033[0m" );
      }

      if ( criteria == 0 )
      {
         return true;
      }
      else
      {
         if ( criteriaContainer != nullptr )
         {
            *criteriaContainer = criteria;
         }
         return false;
      }
   }

   [[nodiscard]] std::string matchTermination( const uint_t terminationCriteria ) const
   {
      using namespace setupMasks;
      std::stringstream ss;

      if ( terminationCriteria & Iterations )
      {
         ss << "IT,";
      }
      if ( terminationCriteria & ResidualThreshold )
      {
         ss << "RT,";
      }
      if ( terminationCriteria & ResidualRate )
      {
         ss << "RR,";
      }
      if ( terminationCriteria & ErrorThreshold )
      {
         ss << "ET,";
      }
      if ( terminationCriteria & ErrorRate )
      {
         ss << "ER,";
      }
      if ( terminationCriteria & ErrIsNaN )
      {
         ss << "\033[33m"
            << "EN,"
            << "\033[0m";
      }
      if ( terminationCriteria & ResIsNaN )
      {
         ss << "\033[33m"
            << "RN,"
            << "\033[0m";
      }

      std::string terminationString = ss.str();
      if ( !terminationString.empty() )
      {
         terminationString.pop_back();
      }
      else
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "No termination criteria was meet but the outer loop terminated anyways." );
      }

      return terminationString;
   }

   const Problem_t problemSetup = {};

}; // struct Setup

template < typename Norm_t >
static constexpr Norm_t chooseRightNorm( const uint_t  normType,
                                         const Norm_t& relativeNorm,
                                         const Norm_t& pointWiseNorm,
                                         const Norm_t& relativeWeakNorm = {} )
{
   using namespace setupMasks;

   if ( normType & Relative || normType & Vector )
   {
      return relativeNorm;
   }
   else if ( normType & Absolute || normType & NormalizedVector )
   {
      return pointWiseNorm;
   }
   else if ( normType & Weak )
   {
      return relativeWeakNorm;
   }
   else
   {
      WALBERLA_ABORT( "No supported norm type was chosen." )
   }
}

} // namespace hyteg::mixedPrecisionMT
