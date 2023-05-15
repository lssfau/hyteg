/*
 * Copyright (c) 2022-2023 Daniel Bauer.
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

#include "hyteg/n1e1functionspace/HybridSmoother.hpp"

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_linear_form_affine_q6.hpp"
#include "hyteg/forms/form_hyteg_manual/N1E1FormCurlCurl.hpp"
#include "hyteg/forms/form_hyteg_manual/N1E1FormMass.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"

#include "common.hpp"

using namespace hyteg;

void test( const uint_t        level,
           const n1e1::System& system,
           const uint_t        nIterations,
           const real_t        expectedResidual,
           const real_t        expectedError,
           const bool          writeVTK = false )
{
   using namespace n1e1;

   using P1LaplaceOperator = P1ConstantLaplaceOperator;
   using N1E1Smoother      = ChebyshevSmoother< N1E1ElementwiseLinearCombinationOperator >;
   using P1Smoother        = GaussSeidelSmoother< P1LaplaceOperator >;

   const uint_t minLevel = level;
   const uint_t maxLevel = level;

   SetupPrimitiveStorage setupStorage( system.domain_, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   writeDomainPartitioningVTK( storage, "../../output", "HybridSmootherTest_partitioning" );

   N1E1VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > f( "f", storage, minLevel, maxLevel );
   N1E1VectorFunction< real_t > sol( "sol", storage, level, level );
   N1E1VectorFunction< real_t > err( "err", storage, level, level );
   N1E1VectorFunction< real_t > tmp( "tmp", storage, level, level );

   // operators
   N1E1Form_curl_curl                       curlCurlForm;
   N1E1Form_mass                            massForm;
   N1E1ElementwiseMassOperator              M( storage, level, level );
   N1E1ElementwiseLinearCombinationOperator A( storage, minLevel, maxLevel, { { 1.0, 1.0 }, { &curlCurlForm, &massForm } } );
   auto p1LaplaceOperator = std::make_shared< P1LaplaceOperator >( storage, minLevel, maxLevel );

   // assemble RHS
   assembleLinearForm< forms::n1e1_linear_form_affine_q6 >( level, level, { system.rhs_ }, f );

   // smoothers
   auto chebyshevSmoother = std::make_shared< N1E1Smoother >( storage, minLevel, maxLevel );

   sol.interpolate( system.analyticalSol_, level );
   const real_t spectralRadius = chebyshev::estimateRadius( A, level, 40, storage, sol, tmp );
   sol.interpolate( system.analyticalSol_, level );

   chebyshevSmoother->setupCoefficients( 4, spectralRadius );
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( spectralRadius );

   std::shared_ptr< Solver< N1E1ElementwiseLinearCombinationOperator > > n1e1Smoother = chebyshevSmoother;
   std::shared_ptr< Solver< P1LaplaceOperator > >                        p1Smoother   = std::make_shared< P1Smoother >();

   // solve Au = f with hybrid smoother
   HybridSmoother hybridSmoother( storage, p1LaplaceOperator, n1e1Smoother, p1Smoother, minLevel, maxLevel );

   A.apply( u, tmp, level, DoFType::Inner );
   tmp.assign( { 1.0, -1.0 }, { f, tmp }, level, DoFType::Inner );
   real_t prevResidual = std::sqrt( tmp.dotGlobal( tmp, level ) );
   // WALBERLA_LOG_DEVEL_VAR_ON_ROOT( prevResidual );

   err.assign( { 1.0, -1.0 }, { sol, u }, level );
   M.apply( err, tmp, level, DoFType::All );
   real_t prevL2Error = std::sqrt( err.dotGlobal( tmp, level ) );
   // WALBERLA_LOG_DEVEL_VAR_ON_ROOT( prevL2Error )

   VTKOutput vtk( "../../output/", "HybridSmootherTest", storage );
   if ( writeVTK )
   {
      vtk.add( u );
      vtk.add( f );
      vtk.add( sol );
      vtk.add( err );
      vtk.write( level, 0 );
   }

   for ( uint_t i = 0; i < nIterations; ++i )
   {
      hybridSmoother.solve( A, u, f, level );

      // determine residual
      A.apply( u, tmp, level, DoFType::Inner );
      tmp.assign( { 1.0, -1.0 }, { f, tmp }, level, DoFType::Inner );
      const real_t residual = std::sqrt( tmp.dotGlobal( tmp, level ) );
      // WALBERLA_LOG_DEVEL_VAR_ON_ROOT( residual )

      WALBERLA_CHECK_LESS( residual, prevResidual, "hybrid smoother, does not decrease residual" );
      prevResidual = residual;

      // determine error
      err.assign( { 1.0, -1.0 }, { sol, u }, level );
      M.apply( err, tmp, level, DoFType::All );
      const real_t L2Error = std::sqrt( err.dotGlobal( tmp, level ) );
      // WALBERLA_LOG_DEVEL_VAR_ON_ROOT( L2Error )

      WALBERLA_CHECK_LESS( L2Error, prevL2Error, "hybrid smoother, does not decrease error" );
      prevL2Error = L2Error;

      if ( writeVTK )
      {
         vtk.write( level, i + 1 );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "2-norm of residual: " << prevResidual );
   WALBERLA_LOG_INFO_ON_ROOT( "Approximate L2-norm of error: " << prevL2Error );

   // the calculated residual should up to 10% match the expected residual
   WALBERLA_CHECK_LESS( 0.9 * expectedResidual, prevResidual, "residual for hybrid smoother has changed " );
   WALBERLA_CHECK_GREATER( 1.1 * expectedResidual, prevResidual, "residual for hybrid smoother has changed " );

   // the calculated error should up to 10% match the expected error
   WALBERLA_CHECK_LESS( 0.9 * expectedError, prevL2Error, "error for hybrid smoother has changed " );
   WALBERLA_CHECK_GREATER( 1.1 * expectedError, prevL2Error, "error for hybrid smoother has changed " );
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "### Single Macro ###" );
   test( 4, n1e1::System::polynomialOnSingleTet(), 60, 0.00317698, 0.000393222 );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "### Pyramid ###" );
   test( 4, n1e1::System::polynomialOnPyramid(), 40, 0.287105, 0.0664059 );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "### Cube ###" );
   test( 3, n1e1::System::sinusoidalOnCube(), 40, 0.59872, 0.060057 );

   return EXIT_SUCCESS;
}
