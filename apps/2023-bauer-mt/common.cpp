/*
* Copyright (c) 2023 Daniel Bauer.
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

#include "common.hpp"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_curl_curl_affine_q0.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_linear_form_affine_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_mass_affine_qe.hpp"
#include "hyteg/gridtransferoperators/N1E1toN1E1Prolongation.hpp"
#include "hyteg/gridtransferoperators/N1E1toN1E1Restriction.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "mixed_operator/HybridSmoother.hpp"

Results solve( const Params& params, const bool useGmg )
{
   using namespace n1e1;

   using P1LaplaceOperator = P1ConstantLinearCombinationOperator;
   using P1LaplaceForm     = P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise >;
   using N1E1Smoother      = ChebyshevSmoother< N1E1ElementwiseLinearCombinationOperator >;
   using P1Smoother        = GaussSeidelSmoother< P1LaplaceOperator >;

   // Mesh
   SetupPrimitiveStorage setupStorage( params.system.domain_, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   if ( params.writeVTK )
   {
      writeDomainPartitioningVTK( storage, "output", params.name + "-partitioning" );
   }

   // Operators
   forms::n1e1_curl_curl_affine_q0 curlCurlForm;
   forms::n1e1_mass_affine_qe      massForm;

   N1E1ElementwiseMassOperator              M( storage, params.minLevel, params.maxLevel );
   N1E1ElementwiseLinearCombinationOperator A(
       storage, params.minLevel, params.maxLevel, { params.coefficients, { &curlCurlForm, &massForm } } );

   if ( params.computeAndStoreLocalElementMatrices )
   {
      A.computeAndStoreLocalElementMatrices();
   }

   // Functions
   N1E1VectorFunction< real_t > u( "u", storage, params.minLevel, params.maxLevel );
   N1E1VectorFunction< real_t > f( "f", storage, params.minLevel, params.maxLevel );
   N1E1VectorFunction< real_t > sol( "sol", storage, params.minLevel, params.maxLevel );
   N1E1VectorFunction< real_t > res( "res", storage, params.minLevel, params.maxLevel );
   N1E1VectorFunction< real_t > err( "err", storage, params.minLevel, params.maxLevel );
   N1E1VectorFunction< real_t > tmp( "tmp", storage, params.minLevel, params.maxLevel );

   // Assemble RHS.
   assembleLinearForm< forms::n1e1_linear_form_affine_q6 >( params.maxLevel, params.maxLevel, { params.system.rhs_ }, f );

   // Initialize
   if ( params.initialGuess.has_value() )
   {
      u.interpolate( params.initialGuess.value(), params.maxLevel, DoFType::Inner );
   }
   // Boundary conditions: homogeneous tangential trace
   u.interpolate( Point3D{ 0.0, 0.0, 0.0 }, params.maxLevel, DoFType::Boundary );

   // Hybrid smoother
   P1LaplaceForm laplaceForm;
   auto          p1LaplaceOperator = std::make_shared< P1ConstantLinearCombinationOperator >(
       storage, params.minLevel, params.maxLevel, P1LinearCombinationForm{ { params.coefficients[1] }, { &laplaceForm } } );

   auto chebyshevSmoother = std::make_shared< N1E1Smoother >( storage, params.minLevel, params.maxLevel );
   auto p1Smoother        = std::make_shared< P1Smoother >();

   // estimate spectral radius (initial guess is sol + u_0)
   sol.interpolate( params.system.analyticalSol_, params.maxLevel );
   sol.assign( { 1.0, 1.0 }, { sol, u }, params.maxLevel );

   storage->getTimingTree()->start( "Chebyshev estimate radius" );
   const real_t spectralRadius =
       chebyshev::estimateRadius( A, params.maxLevel, params.numSpectralRadiusEstIts, storage, sol, tmp );
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( spectralRadius )
   storage->getTimingTree()->stop( "Chebyshev estimate radius" );

   chebyshevSmoother->setupCoefficients(
       params.chebyshevOrder, params.lowerBoundFactor * spectralRadius, params.upperBoundFactor * spectralRadius );

   auto hybridSmoother =
       std::make_shared< HybridSmoother< N1E1ElementwiseLinearCombinationOperator, P1LaplaceOperator > >( storage,
                                                                                                          p1LaplaceOperator,
                                                                                                          chebyshevSmoother,
                                                                                                          p1Smoother,
                                                                                                          params.minLevel,
                                                                                                          params.maxLevel,
                                                                                                          params.n1e1SmoothSteps,
                                                                                                          params.p1SmoothSteps );

   std::shared_ptr< Solver< N1E1ElementwiseLinearCombinationOperator > > solver;
   if ( useGmg )
   {
      // GMG solver
#ifdef HYTEG_BUILD_WITH_PETSC
      WALBERLA_LOG_INFO_ON_ROOT( "Using PETSc solver" )
      auto coarseGridSolver =
          std::make_shared< PETScCGSolver< N1E1ElementwiseLinearCombinationOperator > >( storage, params.minLevel );
#else
      WALBERLA_LOG_INFO_ON_ROOT( "Using HyTeG solver" )
      auto coarseGridSolver = std::make_shared< CGSolver< N1E1ElementwiseLinearCombinationOperator > >(
          storage, params.minLevel, params.minLevel, 10000, 1e-12 );
#endif
      auto restrictionOperator  = std::make_shared< N1E1toN1E1Restriction >();
      auto prolongationOperator = std::make_shared< N1E1toN1E1Prolongation >();

      auto gmgSolver =
          std::make_shared< GeometricMultigridSolver< N1E1ElementwiseLinearCombinationOperator > >( storage,
                                                                                                    hybridSmoother,
                                                                                                    coarseGridSolver,
                                                                                                    restrictionOperator,
                                                                                                    prolongationOperator,
                                                                                                    params.minLevel,
                                                                                                    params.maxLevel,
                                                                                                    params.preSmoothSteps,
                                                                                                    params.postSmoothSteps );
      solver = gmgSolver;
   }
   else
   {
      solver = hybridSmoother;
   }

   // Interpolate solution
   sol.interpolate( params.system.analyticalSol_, params.maxLevel );

   // Determine initial solution norm
   const real_t initU2 = std::sqrt( u.dotGlobal( u, params.maxLevel, DoFType::Inner ) );

   // Determine initial residual
   A.apply( u, tmp, params.maxLevel, DoFType::All );
   res.assign( { 1.0, -1.0 }, { f, tmp }, params.maxLevel, DoFType::Inner );
   const real_t initRes = std::sqrt( res.dotGlobal( res, params.maxLevel, DoFType::Inner ) );

   // Solve system.
   real_t uNorm    = initU2;
   real_t residual = initRes;
   real_t discrL2  = 0.0;
   uint_t its      = 0;

   for ( ;
         its < params.nMaxIterations && ( params.u2Reduction.has_value() ? uNorm / initU2 > params.u2Reduction.value() : true ) &&
         ( params.residual2Reduction.has_value() ? residual / initRes > params.residual2Reduction.value() : true );
         ++its )
   {
      solver->solve( A, u, f, params.maxLevel );

      // Determine solution norm
      uNorm = std::sqrt( u.dotGlobal( u, params.maxLevel, DoFType::Inner ) );

      // Determine residual
      A.apply( u, tmp, params.maxLevel, DoFType::All );
      res.assign( { 1.0, -1.0 }, { f, tmp }, params.maxLevel, DoFType::Inner );
      residual = std::sqrt( res.dotGlobal( res, params.maxLevel, DoFType::Inner ) );
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( residual / initRes )

      // Determine L2 error
      err.assign( { 1.0, -1.0 }, { u, sol }, params.maxLevel );
      M.apply( err, tmp, params.maxLevel, DoFType::All );
      discrL2 = std::sqrt( err.dotGlobal( tmp, params.maxLevel ) );
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT( discrL2 )
   }

   if ( params.writeVTK )
   {
      VTKOutput vtk( "output", params.name, storage );
      vtk.add( u );
      vtk.add( f );
      vtk.add( sol );
      vtk.add( res );
      vtk.add( err );
      vtk.write( params.maxLevel );
   }

   storage->getTimingTree()->synchronize();
   auto ttReduced = storage->getTimingTree()->getReduced().getCopyWithRemainder();

   return Results{
       numberOfGlobalDoFs( u, params.maxLevel ), spectralRadius, initU2, uNorm, initRes, residual, discrL2, its, ttReduced };
}

Results solve( const Params& params )
{
   return solve( params, true );
}

Results smooth( const Params& params )
{
   return solve( params, false );
}
