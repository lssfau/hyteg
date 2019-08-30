
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
        const uint_t numGlobalDoFs = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *u.u.getStorage(), level );

        /////////
        // R_0 //
        /////////

        // velocity
        R.u.assign( {1.0}, {b.u}, level, flag_ );
        R.v.assign( {1.0}, {b.v}, level, flag_ );
        stokesOperator.A.apply( u.u, tmp.u, level, flag_ );
        stokesOperator.A.apply( u.v, tmp.v, level, flag_ );
        stokesOperator.divT_x.apply( u.p, tmp.u, level, flag_, Add );
        stokesOperator.divT_y.apply( u.p, tmp.v, level, flag_, Add );
        R.u.add( {-1.0}, {tmp.u}, level, flag_ );
        R.v.add( {-1.0}, {tmp.v}, level, flag_ );

        // pressure
        stokesOperator.div_x.apply( u.u, tmp.p, level, flag_ | DirichletBoundary, Replace );
        stokesOperator.div_y.apply( u.v, tmp.p, level, flag_ | DirichletBoundary, Add );
        R.p.assign( {-1.0}, {tmp.p}, level, flag_ | DirichletBoundary );

        /////////////
        // R_hat_0 //
        /////////////

        // velocity
        velocityBlockSolver_->solve( stokesOperator.A, RHat.u, R.u, level );
        velocityBlockSolver_->solve( stokesOperator.A, RHat.v, R.v, level );

        // pressure
        stokesOperator.div_x.apply( RHat.u, RHat.p, level, flag_ | DirichletBoundary, Replace );
        stokesOperator.div_y.apply( RHat.v, RHat.p, level, flag_ | DirichletBoundary, Add );
        stokesOperator.div_x.apply( u.u, RHat.p, level, flag_ | DirichletBoundary, Add );
        stokesOperator.div_y.apply( u.v, RHat.p, level, flag_ | DirichletBoundary, Add );

        ///////////////
        // R_tilde_0 //
        ///////////////

        RTilde.assign( {1.0}, {RHat}, level, flag_ );
        invMass.apply( RHat.p, RTilde.p, level, flag_ | DirichletBoundary );

        ///////
        // P //
        ///////

        P.assign( {1.0}, {RTilde}, level, flag_ );

        ///////////
        // alpha //
        ///////////

        alpha_n = RTilde.p.dotGlobal( RHat.p, level, flag_ | DirichletBoundary );

        // calculate M * P
        // 1. calc tmp := A^-1 * B^T * P_pressure
        stokesOperator.divT_x.apply( P.p, tmp.u, level, flag_ );
        stokesOperator.divT_y.apply( P.p, tmp.v, level, flag_ );
        velocityBlockSolver_->solve( stokesOperator.A, tmp_2.u, tmp.u, level );
        velocityBlockSolver_->solve( stokesOperator.A, tmp_2.v, tmp.v, level );
        // 2. velocity of M * P
        MP.u.assign( {1.0, 1.0}, {P.u, tmp_2.u}, level, flag_ );
        MP.v.assign( {1.0, 1.0}, {P.v, tmp_2.v}, level, flag_ );
        // 3. pressure of M * P
        stokesOperator.div_x.apply( tmp_2.u, MP.p, level, flag_ | DirichletBoundary, Replace );
        stokesOperator.div_y.apply( tmp_2.v, MP.p, level, flag_ | DirichletBoundary, Add );

        alpha_d = P.p.dotGlobal( MP.p, level, flag_ | DirichletBoundary );

        alpha = alpha_n / alpha_d;

        // WALBERLA_LOG_INFO_ON_ROOT( "alpha_0 = " << alpha_n << " / " << alpha_d << " = " << alpha );

        ///////
        // X //
        ///////

        // u.add( {alpha}, {P}, level, flag_ );
        u.u.add( {alpha}, {P.u}, level, flag_ );
        u.v.add( {alpha}, {P.v}, level, flag_ );
        u.p.add( {alpha}, {P.p}, level, flag_ | DirichletBoundary );

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
        invMass.apply( RHat.p, RTilde.p, level, flag_ | DirichletBoundary );

        for ( uint_t i = 0; i < maxIterations_; i++ )
        {
            // rescale pressure
            vertexdof::projectMean( u.p, level );

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


            beta_n = RTilde.p.dotGlobal( RHat.p, level, flag_ | DirichletBoundary );

            if ( !solveExactly_ )
            {
              stokesOperator.A.apply( RHat.u, tmp.u, level, flag_ );
              stokesOperator.A.apply( RHat.v, tmp.v, level, flag_ );
              tmp.u.add( {-1.0}, {R.u}, level );
              tmp.v.add( {-1.0}, {R.v}, level );
              beta_n += RHat.u.dotGlobal( tmp.u, level, flag_ );
              beta_n += RHat.v.dotGlobal( tmp.v, level, flag_ );
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
            stokesOperator.divT_x.apply( P.p, tmp.u, level, flag_ );
            stokesOperator.divT_y.apply( P.p, tmp.v, level, flag_ );
            velocityBlockSolver_->solve( stokesOperator.A, tmp_2.u, tmp.u, level );
            velocityBlockSolver_->solve( stokesOperator.A, tmp_2.v, tmp.v, level );
            // 2. velocity of M * P
            MP.u.assign( {1.0, 1.0}, {P.u, tmp_2.u}, level, flag_ );
            MP.v.assign( {1.0, 1.0}, {P.v, tmp_2.v}, level, flag_ );
            // 3. pressure of M * P
            stokesOperator.div_x.apply( tmp_2.u, MP.p, level, flag_ | DirichletBoundary, Replace );
            stokesOperator.div_y.apply( tmp_2.v, MP.p, level, flag_ | DirichletBoundary, Add );

            alpha_d = P.p.dotGlobal( MP.p, level, flag_ | DirichletBoundary );

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

            u.u.add( {alpha}, {P.u}, level, flag_ );
            u.v.add( {alpha}, {P.v}, level, flag_ );
            u.p.add( {alpha}, {P.p}, level, flag_ | DirichletBoundary );

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
            invMass.apply( RHat.p, RTilde.p, level, flag_ | DirichletBoundary );

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

