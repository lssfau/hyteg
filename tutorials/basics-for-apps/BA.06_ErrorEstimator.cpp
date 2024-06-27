/*
 * Copyright (c) 2024 Benjamin Mann, Marcus Mohr.
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

#include "hyteg/adaptiverefinement/error_estimator.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::math::pi;

namespace hyteg {

/**
 * \page BA.06_ErrorEstimator Tutorial BA.06 - Using an ErrorEstimator for AMR
 *
 * \dontinclude tutorials/basics-for-apps/BA.06_ErrorEstimator.cpp
 *
 * \brief In this tutorial we will go through all required steps to set up an ErrorEstimator and apply it for adaptive refinement.
 *
 * \section T06-intro Introduction
 *
 * When solving a problem, we might want to have an idea whether the resolution of our mesh is sufficiently fine.
 * Since, in practice, we don't have access to the actual error, we require an error estimator.
 * This is particularly true if the problem demands adaptive refinement.
 * In this case, we not only need an estimate to the global error, but also some indication of where in the domain the errors are largest.
 *
 * The class ErrorEstimator provides such functionality.
 * It is intended to be used with adaptiveRefinement::Mesh, but can also be applied for non adaptive setups.
 *
 * In the following example we set up a full app solving a problem Au=b.
 * We assume that the problem requires a heterogeneous resolution and therefore iteratively apply AMR.
 * In each iteration, we solve the problem with a more adapted macro grid.
 *
 * \section T06-adaptiveMesh Creating an adaptive mesh
 *
 * First, we follow the steps introduced in \ref BA.05_AMR to create an adaptiveRefinement::Mesh.
 *
 * \snippet{trimleft} this AdaptiveMesh
 *
 * We then run a loop over the following steps.
 *
 * \section T06-problem Problem definition on the current grid
 *
 * We first define the problem as described in previous tutorials
 *
 * \snippet{trimleft} this Problem
 *
 * \section T06-errest Error estimate
 * Our ErrorEstimator is defined as follows:
 *
 * Let \f$ \ell=(L-j) \f$ for some \f$ j \in \{1, ...,  (L-\ell_{min}-1) \} \f$. Then the discretization error on level \f$\ell\f$ can be estimated by
 * \f[
   \begin{equation*}
   \tilde{e}_\ell = u_L - u_\ell
   \end{equation*}
 * \f]
 * Ultimately, we are interested in the error at level L, though.
 * Therefore, we use the lower level errors to estimate the convergence factor by
 * \f[
   \begin{equation*}
   \frac{\|\tilde{e}_\ell\|}{\|\tilde{e}_{\ell-1}\|}
   \end{equation*}
 * \f]
 *
 * We then scale \f$\|\tilde{e}_\ell\|\f$ by the appropriate power of this estimated convergence factor to yield an estimate to \f$\|e_L\|\f$.
 * Therefore, the parameter \f$ j\f$ influences the accuracy and robustness of the obtained estimate:
 * Larger \f$ j\f$ yield a more accurate estimate, while smaller \f$ j\f$ makes it less sensitive to inaccuracies in the convergence estimate due to pre-asymptotic convergence.
 * In practice \f$ j=2\f$ seems to be a good compromise.
 *
 * The constructor of our ErrorEstimator takes an argument `j_max`.
 * We can later evaluate the above defined estimator for all \f$ j \leq j_\text{max}\f$.
 *
 * \snippet{trimleft} this ErrorEstimator
 *
 * \section T06-solver FMG Solver
 * Next, we define the solver.
 * Due to the above definition, we need solutions \f$ u_\ell \f$ from coarser levels for our estimate.
 * Therefore, a full multigrid solver must be used.
 * It is also required, that the function ErrorEstimator::fmg_callback() is given to the constructor of FullMultigridSolver as postProlongateCallback.
 *
 * \snippet{trimleft} this Solver
 *
 * \section T06-solve Solving the problem and estimating the global error
 *
 * After solving the problem using FMG, we call `ErrorEstimator::estimate` to compute the estimates for all \f$ j=1,...,j_{max}\f$.
 *
 * \snippet{trimleft} this Estimate
 *
 * Now, we break the refinement loop if the estimated error is lower than some threshold
 *
 * \snippet{trimleft} this GlobalError
 *
 * \section T06-refine Refining the mesh based on the local error indicator
 *
 * For the adaptive refinement, we require an indication of the local errors.
 * The method `ErrorEstimator::eta_T_sq` returns a vector of estimated squared local errors \f$ \|\tilde{e}_{L-1}\|_{L^2(T)}^2 \f$ for all macro faces T (or cells in 3d), giving a good indication where to refine.
 *
 * \snippet{trimleft} this LocalError
 *
 * We now must choose a refinement strategy.
 * Currently, there are two strategies available, weighted mean and percentile.
 * The former refines all macro elements where the squared error is greater than the weighted mean squared local error over all elements.
 * The latter refines a given fraction of the macros.
 * Here, we use `RefinementStrategy::WEIGHTED_MEAN`.
 * Setting the weighting parameter to zero yields the arithmetic mean.
 *
 * \snippet{trimleft} this RefinementStrategy
 *
 * Now, we call Mesh::refineRG with the above parameters and provide the local error vector.
 *
 * \snippet{trimleft} this RefineRG
 *
 * The refined mesh is distributed over the processes such that children are owned by the same process as their parents.
 * As this is in general not desirable, `loadbalancing` must be called again.
 * Since we want to continue with the refined mesh in the next iteration, we call `Mesh::make_storage` again, to update the PrimitiveStorage.
 *
 * \snippet{trimleft} this UpdateStorage
 *
 * Here, we overwrite the smart pointer, throwing away the previous storage.
 * Of course, you may also keep it, along with corresponding function data.
 * This may be particularly interesting when planing to interpolate functions from the old to the new grid.
 * Note that doing so requires calling both `make_storage` and the interpolation before `loadbalancing`.
 * This may lead to severe load imbalance during the interpolation step.
 * Therefore, it is currently not recommended to interpolate data between different PrimitiveStorages.
 *
 * \section code Complete Program
 * \include tutorials/basics-for-apps/BA.05_AMR.cpp
 *
 *
 */

void errorEstimatorTutorial()
{
   /// [AdaptiveMesh]
   auto                  meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   adaptiveRefinement::Mesh adaptive_mesh( setupStorage );
   adaptive_mesh.loadbalancing();
   auto storage = adaptive_mesh.make_storage();
   /// [AdaptiveMesh]

   // define min and max lvl for multigrid solver
   uint_t l_min = 0;
   uint_t L     = 5;

   for ( uint_t refinement_step = 0; refinement_step < 100; ++refinement_step )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "REFINEMENT STEP: " << refinement_step );
      WALBERLA_LOG_INFO_ON_ROOT( "number of macro faces: " << adaptive_mesh.n_elements() );

      /// [Problem]
      // problem definition
      using Laplace = P1ConstantLaplaceOperator;
      using P1      = P1Function< real_t >;

      Laplace A( storage, l_min, L );
      P1      u_h( "u_h", storage, l_min, L );
      P1      b( "b", storage, l_min, L );
      // Initialize dirichlet data
      auto g = []( const hyteg::Point3D& x ) { return std::sin( 32 * pi * x[0] * x[0] ) * std::cos( 32 * pi * x[1] * x[1] ); };
      for ( auto lvl = l_min; lvl <= L; ++lvl )
         u_h.interpolate( g, lvl, DirichletBoundary );
      /// [Problem]

      /// [ErrorEstimator]
      // define the errorestimator
      uint_t                                   j_max = 2;
      adaptiveRefinement::ErrorEstimator< P1 > errorEstimator( u_h, j_max );
      /// [ErrorEstimator]

      /// [Solver]
      auto P   = std::make_shared< P1toP1LinearProlongation<> >();
      auto R   = std::make_shared< P1toP1LinearRestriction<> >();
      auto S   = std::make_shared< GaussSeidelSmoother< Laplace > >();
      auto cg  = std::make_shared< CGSolver< Laplace > >( storage, l_min, l_min );
      auto gmg = std::make_shared< GeometricMultigridSolver< Laplace > >( storage, S, cg, R, P, l_min, L );

      FullMultigridSolver< Laplace > fmg( storage, gmg, P, l_min, L, 1, []( uint_t ) {}, errorEstimator.fmg_callback() );
      /// [Solver]

      /// [Estimate]
      // solve the problem
      fmg.solve( A, u_h, b, L );
      // compute local and global error estimates
      errorEstimator.estimate();
      /// [Estimate]

      /// [GlobalError]
      // estimate global L2 error
      uint_t j   = 2; //  1 <= j <= j_max
      real_t eta = errorEstimator.eta( j );
      WALBERLA_LOG_INFO_ON_ROOT( "||e||_L2 ≈ eta_j = " << eta );
      // stopping criterion
      if ( eta < 1e-4 )
         break;
      /// [GlobalError]

      /// [LocalError]
      // get an indication to local squared L2 errors
      auto local_sqared_err = errorEstimator.eta_T_sq();
      /// [LocalError]

      /// [RefinementStrategy]
      // choose a refinement strategy
      auto   strategy  = adaptiveRefinement::RefinementStrategy::WEIGHTED_MEAN;
      real_t weighting = 0;
      bool   verbose   = true;
      /// [RefinementStrategy]

      /// [RefineRG]
      // refine those macros where the squared local errors are largest
      adaptive_mesh.refineRG( local_sqared_err, strategy, weighting, verbose );
      /// [RefineRG]

      /// [UpdateStorage]
      // get storage from refined mesh for next iteration
      adaptive_mesh.loadbalancing();
      storage = adaptive_mesh.make_storage();
      /// [UpdateStorage]
   }
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   hyteg::errorEstimatorTutorial();
   return 0;
}
