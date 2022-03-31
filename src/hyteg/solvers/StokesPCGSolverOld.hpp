/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#pragma once

#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "core/debug/CheckFunctions.h"

namespace hyteg {

///  \brief Preconditioned CG solver for the Stokes system
///         compare w/ Elman (1996): Multigrid and Krylov subspace methods, sect. 2.2
template< class OperatorType >
class StokesPCGSolverOld : public Solver< OperatorType >
{
public:

    typedef typename OperatorType::srcType FunctionType;  

    StokesPCGSolverOld( const std::shared_ptr< PrimitiveStorage > & storage,
                   const std::shared_ptr< Solver< typename OperatorType::VelocityOperator_T > > & velocityBlockSolver,
                   const uint_t & minLevel, 
                   const uint_t & maxLevel,
                   const real_t & residualTolerance, 
                   const uint_t & maxIterations, 
                   const DoFType & flag_,
                   const bool & clipAlpha = true,
                   const bool & clipBeta  = false ) :
            velocityBlockSolver_( velocityBlockSolver ),
            minLevel_( minLevel ),
            maxLevel_( maxLevel ),
            residualTolerance_( residualTolerance ),
            maxIterations_( maxIterations ),
            flag_( flag_ ),

            clipAlpha_( clipAlpha ),
            clipBeta_( clipBeta ),

            tmp("tmp", storage, minLevel, maxLevel),
            tmp_2("tmp_2", storage, minLevel, maxLevel),
            tmp_3("tmp_3", storage, minLevel, maxLevel),
            z("z", storage, minLevel, maxLevel),
            s("s", storage, minLevel, maxLevel),
            v("v", storage, minLevel, maxLevel),
            BTs("BTs", storage, minLevel, maxLevel),
            Bv("Bv", storage, minLevel, maxLevel),

            R("R", storage, minLevel, maxLevel),
            RHat("RHat", storage, minLevel, maxLevel),
            RTilde("RTilde", storage, minLevel, maxLevel),
            P("P", storage, minLevel, maxLevel),
            MP("SchurP", storage, minLevel, maxLevel),
            Au("Au", storage, minLevel, maxLevel),
            residual("residualP", storage, minLevel, maxLevel),

            invMass( storage, minLevel, maxLevel )

    {
        const bool available = !storage->hasGlobalCells() && std::is_same< OperatorType, P2P1TaylorHoodStokesOperator >::value;
        WALBERLA_CHECK( available, "SchurCG solver currently only available for Taylor-Hood discretization on 2D domains." );
    }

    void solve( const OperatorType & stokesOperator,
                const FunctionType & u,
                const FunctionType & b,
                const uint_t level )
    {
        real_t alpha, alpha_n, alpha_d, beta, beta_n, beta_d;
        const uint_t numGlobalDoFs = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *u.uvw()[0].getStorage(), level );

        /////////
        // R_0 //
        /////////

        // velocity
        R.uvw()[0].assign( {1.0}, {b.uvw()[0]}, level, flag_ );
        R.uvw()[1].assign( {1.0}, {b.uvw()[1]}, level, flag_ );
        stokesOperator.getA().apply( u.uvw()[0], tmp.uvw()[0], level, flag_ );
        stokesOperator.getA().apply( u.uvw()[1], tmp.uvw()[1], level, flag_ );
        stokesOperator.divT.apply( u.p(), tmp.uvw(), level, flag_, Add );
        R.uvw()[0].add( {-1.0}, {tmp.uvw()[0]}, level, flag_ );
        R.uvw()[1].add( {-1.0}, {tmp.uvw()[1]}, level, flag_ );

        // pressure
        stokesOperator.div.apply( u.uvw(), tmp.p(), level, flag_ | DirichletBoundary, Replace );
        R.p().assign( {-1.0}, {tmp.p()}, level, flag_ | DirichletBoundary );

        /////////////
        // R_hat_0 //
        /////////////

        // velocity
        velocityBlockSolver_->solve( stokesOperator.getA(), RHat.uvw()[0], R.uvw()[0], level );
        velocityBlockSolver_->solve( stokesOperator.getA(), RHat.uvw()[1], R.uvw()[1], level );

        // pressure
        stokesOperator.div.apply( RHat.uvw(), RHat.p(), level, flag_ | DirichletBoundary, Replace );
        stokesOperator.div.apply( u.uvw(), RHat.p(), level, flag_ | DirichletBoundary, Add );

        ///////////////
        // R_tilde_0 //
        ///////////////

        RTilde.assign( {1.0}, {RHat}, level, flag_ );
        invMass.apply( RHat.p(), RTilde.p(), level, flag_ | DirichletBoundary );

        ///////
        // P //
        ///////

        P.assign( {1.0}, {RTilde}, level, flag_ );

        ///////////
        // alpha //
        ///////////

        alpha_n = RTilde.p().dotGlobal( RHat.p(), level, flag_ | DirichletBoundary );

        // calculate M * P
        // 1. calc tmp := A^-1 * B^T * P_pressure
        stokesOperator.divT.apply( P.p(), tmp.uvw(), level, flag_ );
        velocityBlockSolver_->solve( stokesOperator.getA(), tmp_2.uvw()[0], tmp.uvw()[0], level );
        velocityBlockSolver_->solve( stokesOperator.getA(), tmp_2.uvw()[1], tmp.uvw()[1], level );
        // 2. velocity of M * P
        MP.uvw()[0].assign( {1.0, 1.0}, {P.uvw()[0], tmp_2.uvw()[0]}, level, flag_ );
        MP.uvw()[1].assign( {1.0, 1.0}, {P.uvw()[1], tmp_2.uvw()[1]}, level, flag_ );
        // 3. pressure of M * P
        stokesOperator.div.apply( tmp_2.uvw(), MP.p(), level, flag_ | DirichletBoundary, Replace );

        alpha_d = P.p().dotGlobal( MP.p(), level, flag_ | DirichletBoundary );

        alpha = alpha_n / alpha_d;

        // WALBERLA_LOG_INFO_ON_ROOT( "alpha_0 = " << alpha_n << " / " << alpha_d << " = " << alpha );

        ///////
        // X //
        ///////

        // u.add( {alpha}, {P}, level, flag_ );
        u.uvw()[0].add( {alpha}, {P.uvw()[0]}, level, flag_ );
        u.uvw()[1].add( {alpha}, {P.uvw()[1]}, level, flag_ );
        u.p().add( {alpha}, {P.p()}, level, flag_ | DirichletBoundary );

        //////////////
        // R update //
        //////////////

        stokesOperator.apply( P, tmp, level, flag_ );
        R.add( {-alpha}, {tmp}, level, flag_ );

        //////////////////
        // R_hat update //
        //////////////////

        RHat.add( {-alpha}, {MP}, level, flag_ );

        ////////////////////
        // R_tilde update //
        ////////////////////

        RTilde.assign( {1.0}, {RHat}, level, flag_ );
        invMass.apply( RHat.p(), RTilde.p(), level, flag_ | DirichletBoundary );

        for ( uint_t i = 0; i < maxIterations_; i++ )
        {
            // rescale pressure
            vertexdof::projectMean( u.p(), level );

            // check residual
            stokesOperator.apply( u, Au, level, flag_ );
            residual.assign( {1.0, -1.0}, {b, Au}, level, flag_ );
            real_t currentResidual = std::sqrt( residual.dotGlobal( residual, level, flag_ ) ) / numGlobalDoFs;
            WALBERLA_LOG_DEVEL_ON_ROOT( "CURRENT SCHUR CG RESIDUAL (" << i << "): " << std::scientific << currentResidual  );
            if ( currentResidual < residualTolerance_ )
              return;

            //////////
            // beta //
            //////////


            beta_n = RTilde.p().dotGlobal( RHat.p(), level, flag_ | DirichletBoundary );

            if ( !solveExactly_ )
            {
              stokesOperator.getA().apply( RHat.uvw()[0], tmp.uvw()[0], level, flag_ );
              stokesOperator.getA().apply( RHat.uvw()[1], tmp.uvw()[1], level, flag_ );
              tmp.uvw()[0].add( {-1.0}, {R.uvw()[0]}, level );
              tmp.uvw()[1].add( {-1.0}, {R.uvw()[1]}, level );
              beta_n += RHat.uvw()[0].dotGlobal( tmp.uvw()[0], level, flag_ );
              beta_n += RHat.uvw()[1].dotGlobal( tmp.uvw()[1], level, flag_ );
            }

            beta_d = alpha_n;

            if ( clipBeta_ )
            {
              if ( beta_n < 1e-16 && beta_d < 1e-16 )
              {
                beta = 1.0;
              }
              else
              {
                beta = beta_n / beta_d;
                beta = beta > 1.0 ? 1.0 : beta;
              }
            }
            else
            {
              beta = beta_n / beta_d;
            }
            // WALBERLA_LOG_INFO_ON_ROOT( "beta_" << i << " = " << beta_n << " / " << beta_d << " = " << beta );

            //////////////
            // P update //
            //////////////

            P.assign( {1.0, beta}, {RTilde, P}, level, flag_ );

            ///////////
            // alpha //
            ///////////

            alpha_n = beta_n;

            // calculate M * P
            // 1. calc tmp := A^-1 * B^T * P_pressure
            stokesOperator.divT.apply( P.p(), tmp.uvw(), level, flag_ );
            velocityBlockSolver_->solve( stokesOperator.getA(), tmp_2.uvw()[0], tmp.uvw()[0], level );
            velocityBlockSolver_->solve( stokesOperator.getA(), tmp_2.uvw()[1], tmp.uvw()[1], level );
            // 2. velocity of M * P
            MP.uvw()[0].assign( {1.0, 1.0}, {P.uvw()[0], tmp_2.uvw()[0]}, level, flag_ );
            MP.uvw()[1].assign( {1.0, 1.0}, {P.uvw()[1], tmp_2.uvw()[1]}, level, flag_ );
            // 3. pressure of M * P
            stokesOperator.div.apply( tmp_2.uvw(), MP.p(), level, flag_ | DirichletBoundary, Replace );

            alpha_d = P.p().dotGlobal( MP.p(), level, flag_ | DirichletBoundary );

            if ( !solveExactly_ )
            {

            }

            if ( clipAlpha_ )
            {
              if ( alpha_n < 1e-16 && alpha_d < 1e-16 )
              {
                alpha = 1.0;
              }
              else
              {
                alpha = alpha_n / alpha_d;
                alpha = alpha > 1.0 ? 1.0 : alpha;
              }
            }
            else
            {
              alpha = alpha_n / alpha_d;
            }
            // WALBERLA_LOG_INFO_ON_ROOT( "alpha_" << i+1 << " = " << alpha_n << " / " << alpha_d << " = " << alpha );

            ///////
            // X //
            ///////

            u.uvw()[0].add( {alpha}, {P.uvw()[0]}, level, flag_ );
            u.uvw()[1].add( {alpha}, {P.uvw()[1]}, level, flag_ );
            u.p().add( {alpha}, {P.p()}, level, flag_ | DirichletBoundary );

            //////////////
            // R update //
            //////////////

            stokesOperator.apply( P, tmp, level, flag_ );
            R.add( {-alpha}, {tmp}, level, flag_ );

            //////////////////
            // R_hat update //
            //////////////////

            RHat.add( {-alpha}, {MP}, level, flag_ );

            ////////////////////
            // R_tilde update //
            ////////////////////

            RTilde.assign( {1.0}, {RHat}, level, flag_ );
            invMass.apply( RHat.p(), RTilde.p(), level, flag_ | DirichletBoundary );

        }

    }

private:
    std::shared_ptr< Solver< typename OperatorType::VelocityOperator_T > > velocityBlockSolver_;
    uint_t minLevel_;
    uint_t maxLevel_;
    real_t residualTolerance_;
    uint_t maxIterations_;
    DoFType flag_;
    bool clipAlpha_;
    bool clipBeta_;
    FunctionType tmp, tmp_2, tmp_3, z, s, v, BTs, Bv,
            R, RHat, RTilde, P, MP, Au, residual;
    P1LumpedInvMassOperator invMass;
    bool solveExactly_;
};

}

