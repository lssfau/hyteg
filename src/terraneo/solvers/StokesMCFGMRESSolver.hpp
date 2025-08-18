/*
 * Copyright (c) 2025 Ponsuganth Ilangovan P
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

#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"

#include "terraneo/helpers/TerraNeoParameters.hpp"
#include "terraneo/solvers/MCSolverBase.hpp"

namespace hyteg {
template < typename StokesOperator_T >
class StokesMCFGMRESSolver : public MCSolverBase< StokesOperator_T >
{
 public:
   // clang-format off
   enum class Parameters
   {
    POWER_ITERATIONS,
        
        FGMRES_OUTER_ITER,
        FGMRES_OUTER_TOL,
            ABLOCK_MG_PRESMOOTH,    
            ABLOCK_MG_POSTSMOOTH,
                ABLOCK_COARSE_ITER,
                ABLOCK_COARSE_TOL,
                ABLOCK_COARSE_GRID_PETSC,
            
            SCHUR_CG_SOLVER_ITER,
            SCHUR_CG_SOLVER_TOL,
   };
   // clang-format on

   StokesMCFGMRESSolver( const std::shared_ptr< PrimitiveStorage >&                 storage,
                         const uint_t                                               minLevel,
                         const uint_t                                               maxLevel,
                         const std::shared_ptr< StokesOperator_T >&                 stokesOperator,
                         const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp1,
                         const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp2,
                         const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp3,
                         terraneo::TerraNeoParameters&                              TN )
   : MCSolverBase< StokesOperator_T >( storage, minLevel, maxLevel )
   , TN_( TN )
   , stokesOperator_( stokesOperator )
   , temp1_( temp1 )
   , temp2_( temp2 )
   , temp3_( temp3 )
   {
      chebyshevABlockSmoother_ =
          std::make_shared< ChebyshevSmoother< typename StokesOperator_T::ViscousOperatorFS_T > >( storage, minLevel, maxLevel );

      real_t spectralRadiusA = 4.0;
      {
         std::function< real_t( const Point3D& ) > randFuncA = []( const Point3D& ) {
            return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
         };

         WALBERLA_LOG_INFO_ON_ROOT( "Estimate spectral radius!" );

         // avoid that the startpoint of our poweriteration is in the kernel of the operator
         temp2->uvw().interpolate( randFuncA, maxLevel, All );
         temp3->uvw().interpolate( randFuncA, maxLevel, All );

         spectralRadiusA = chebyshev::estimateRadius(
             stokesOperator_->getA(), maxLevel, TN.solverParameters.numPowerIterations, storage, temp2->uvw(), temp3->uvw() );

         temp2->uvw().interpolate( 0, maxLevel, All );
         temp3->interpolate( 0, maxLevel, All );

         WALBERLA_LOG_INFO_ON_ROOT( "Estimated spectral radius: " << spectralRadiusA );
      }

      chebyshevABlockSmoother_->setupCoefficients( TN.solverParameters.ChebyshevOrder,
                                                   spectralRadiusA,
                                                   TN.solverParameters.ChebyshevSpectralRadiusUpperLimit,
                                                   TN.solverParameters.ChebyshevSpectralRadiusLowerLimit );
      //   chebyshevABlockSmoother_->setupCoefficients(3u, 4.0, 1.2, 0.1);

      minresCoarseGridABlockSolver_ = std::make_shared< MinResSolver< typename StokesOperator_T::ViscousOperatorFS_T > >(
          storage,
          minLevel,
          maxLevel,
          TN.solverParameters.ABlockCoarseGridIterations,
          TN.solverParameters.ABlockCoarseGridTolerance );

      prolongationOperator_ = std::make_shared< P2toP2QuadraticVectorProlongation >();
      restrictionOperator_  = std::make_shared< P2toP2QuadraticVectorRestriction >();

      gmgABlockSolver_ = std::make_shared< GeometricMultigridSolver< typename StokesOperator_T::ViscousOperatorFS_T > >(
          storage,
          chebyshevABlockSmoother_,
          minresCoarseGridABlockSolver_,
          restrictionOperator_,
          prolongationOperator_,
          minLevel,
          maxLevel,
          TN.solverParameters.ABlockMGPreSmooth,
          TN.solverParameters.ABlockMGPostSmooth );

      cgSchurSolver_ = std::make_shared< CGSolver< typename StokesOperator_T::SchurOperator_T > >(
          storage,
          minLevel,
          maxLevel,
          TN.solverParameters.SchurCoarseGridIterations,
          TN.solverParameters.SchurCoarseGridTolerance );

      blockPreconditioner_ = std::make_shared< BlockFactorisationPreconditioner< StokesOperator_T,
                                                                                 typename StokesOperator_T::ViscousOperatorFS_T,
                                                                                 typename StokesOperator_T::SchurOperator_T > >(
          storage, minLevel, maxLevel, stokesOperator->getSchur(), gmgABlockSolver_, cgSchurSolver_ );

      fgmresFinalSolver_ = std::make_shared< FGMRESSolver< StokesOperator_T > >( storage,
                                                                                 minLevel,
                                                                                 maxLevel,
                                                                                 TN.solverParameters.FGMRESOuterIterations,
                                                                                 TN.solverParameters.FGMRESTolerance,
                                                                                 TN.solverParameters.FGMRESTolerance,
                                                                                 blockPreconditioner_,
                                                                                 TN.solverParameters.FGMRESRestartLength,
                                                                                 TN.solverParameters.FGMRESTolerance );

      fgmresFinalSolver_->setPrintInfo( true );

      this->solverPtr_ = fgmresFinalSolver_;
   }

 private:
   terraneo::TerraNeoParameters& TN_;

   const std::shared_ptr< StokesOperator_T >& stokesOperator_;

   const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp1_;
   const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp2_;
   const std::shared_ptr< P2P1TaylorHoodFunction< real_t > >& temp3_;

   std::shared_ptr< ChebyshevSmoother< typename StokesOperator_T::ViscousOperatorFS_T > > chebyshevABlockSmoother_;
   std::shared_ptr< MinResSolver< typename StokesOperator_T::ViscousOperatorFS_T > >      minresCoarseGridABlockSolver_;

   std::shared_ptr< P2toP2QuadraticVectorProlongation > prolongationOperator_;
   std::shared_ptr< P2toP2QuadraticVectorRestriction >  restrictionOperator_;

   std::shared_ptr< GeometricMultigridSolver< typename StokesOperator_T::ViscousOperatorFS_T > > gmgABlockSolver_;

   std::shared_ptr< CGSolver< typename StokesOperator_T::SchurOperator_T > > cgSchurSolver_;

   std::shared_ptr< BlockFactorisationPreconditioner< StokesOperator_T,
                                                      typename StokesOperator_T::ViscousOperatorFS_T,
                                                      typename StokesOperator_T::SchurOperator_T > >
       blockPreconditioner_;

   std::shared_ptr< FGMRESSolver< StokesOperator_T > > fgmresFinalSolver_;
};
} // namespace hyteg
