
#pragma once

#include "tinyhhg_core/solvers/Solver.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "core/debug/CheckFunctions.h"

namespace hyteg {

///  \brief Preconditioned CG solver for the Stokes system
///         compare w/ Peters 2005, Fast iterative solvers for discrete Stokes equations (sec. 3.1)
template< class OperatorType >
class StokesPCGSolver : public Solver< OperatorType >
{
public:

    typedef typename OperatorType::srcType FunctionType;

    StokesPCGSolver( const std::shared_ptr< PrimitiveStorage > & storage,
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
    z("z", storage, minLevel, maxLevel),
    s("s", storage, minLevel, maxLevel),
    d("d", storage, minLevel, maxLevel),
    q("q", storage, minLevel, maxLevel),
    w("w", storage, minLevel, maxLevel),
    R("R", storage, minLevel, maxLevel),
    RHat("RHat", storage, minLevel, maxLevel),
    P("P", storage, minLevel, maxLevel),
    Au("Au", storage, minLevel, maxLevel),
    residual("residual", storage, minLevel, maxLevel),

    invMass( storage, minLevel, maxLevel ),
    solveExactly_( false )

    {
      const bool available = !storage->hasGlobalCells() && std::is_same< OperatorType, P2P1TaylorHoodStokesOperator >::value;
      WALBERLA_CHECK( available, "SchurCG solver currently only available for Taylor-Hood discretization on 2D domains." );
    }

    void solve( const OperatorType & stokesOperator,
                const FunctionType & u,
                const FunctionType & b,
                const uint_t level )
    {
      real_t alpha, alpha_n, alpha_d, beta, beta_n, beta_d, beta_n_old, beta_n_new, prod0, prod1, prod2;
      prod0 = 0.0;
      prod1 = 0.0;

      const uint_t numGlobalDoFs = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *u.u.getStorage(), level );

      // rHat //
      stokesOperator.apply( u, tmp, level, flag_ );
      RHat.assign( {1.0, -1.0}, {b, tmp}, level );

      // r //

      // apply G to rHat
      velocityBlockSolver_->solve( stokesOperator.A, R.u, RHat.u, level );
      velocityBlockSolver_->solve( stokesOperator.A, R.v, RHat.v, level );
      stokesOperator.div_x.apply( R.u, R.p, level, flag_, Replace );
      stokesOperator.div_y.apply( R.v, R.p, level, flag_, Add );
      R.p.add( {-1.0}, {RHat.p}, level, flag_ );

      for ( uint_t iteration = 0; iteration < maxIterations_; iteration++ )
      {
        // rescale pressure
        vertexdof::projectMean( u.p, level );

        // check residual
        stokesOperator.apply( u, Au, level, flag_ );
        residual.assign( {1.0, -1.0}, {b, Au}, level, flag_ );
        real_t currentResidual = std::sqrt( residual.dotGlobal( residual, level, flag_ ) ) / numGlobalDoFs;
        if ( currentResidual < residualTolerance_ )
          return;

        // z // (preconditioning)
        z.u.assign( {1.0}, {R.u}, level, flag_ );
        z.v.assign( {1.0}, {R.v}, level, flag_ );
        invMass.apply( R.p, z.p, level, flag_ );

        // d //
        stokesOperator.A.apply( R.u, d.u, level, flag_ );
        stokesOperator.A.apply( R.v, d.v, level, flag_ );

        // beta //

        prod0 = d.u.dotGlobal( R.u, level, flag_ ) + d.v.dotGlobal( R.v, level, flag_ );
        prod1 = RHat.u.dotGlobal( R.u, level, flag_ ) + RHat.v.dotGlobal( R.v, level, flag_ );
        prod2 = z.p.dotGlobal( R.p, level, flag_ );
        beta_n_old = beta_n_new;
        beta_n_new = prod0 - prod1 + prod2;
        if ( iteration == 0 )
          beta = 0.0;
        else
          beta = beta_n_new / beta_n_old;

        // p
        P.assign( {1.0, beta}, {z, P}, level, flag_ );

        // s
        s.assign( { 1.0, beta }, { d, s }, level, flag_ );

        // q1 and q2
        // TODO: here C should probably enter q2 (q.p)
        q.u.assign({1.0}, {s.u}, level, flag_);
        q.v.assign({1.0}, {s.v}, level, flag_);
        stokesOperator.divT_x.apply( P.p, q.u, level, flag_, Add );
        stokesOperator.divT_y.apply( P.p, q.v, level, flag_, Add );
        stokesOperator.div_x.apply( P.u, q.p, level, flag_, Replace );
        stokesOperator.div_y.apply( P.v, q.p, level, flag_, Add );

        // w
        velocityBlockSolver_->solve( stokesOperator.A, w.u, q.u, level );
        velocityBlockSolver_->solve( stokesOperator.A, w.v, q.v, level );
        stokesOperator.div_x.apply( w.u, w.p, level, flag_, Replace );
        stokesOperator.div_y.apply( w.v, w.p, level, flag_, Add );
        w.p.add( {-1.0}, {q.p}, level, flag_ );

        // alpha

        prod0 = w.u.dotGlobal( s.u, level, flag_ ) + w.v.dotGlobal( s.v, level, flag_ );
        prod1 = q.u.dotGlobal( P.u, level, flag_ ) + q.v.dotGlobal( P.v, level, flag_ );
        prod2 = w.p.dotGlobal( P.p, level, flag_ );
        alpha_d = prod0 - prod1 + prod2;
        alpha_n = beta_n_new;
        alpha = alpha_n / alpha_d;


        if ( clipAlpha_ )
        {
          if ( alpha_n < 1e-10 && alpha_d < 1e-10 )
          {
            WALBERLA_LOG_DEVEL_ON_ROOT( "== clipping alpha (tolerance)" )
            alpha = 1.0;
          }
          else if ( alpha > 1.0 )
          {
            WALBERLA_LOG_DEVEL_ON_ROOT( "== clipping alpha (greater 1.0)" )
            alpha = 1.0;
          }
        }

        WALBERLA_LOG_DEVEL_ON_ROOT( "CURRENT SCHUR CG RESIDUAL (" << iteration << "): " << std::scientific << currentResidual <<  " | alpha: " << alpha << " | beta: " << beta  );


        // update solution
        u.add( {alpha}, {P}, level, flag_ );

        // r
        R.add( {-alpha}, {w}, level, flag_ );

        // rhat
        RHat.u.add( {-alpha}, {q.u}, level, flag_ );
        RHat.v.add( {-alpha}, {q.v}, level, flag_ );
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
    FunctionType tmp, z, s, d, q, w, R, RHat, P, Au, residual;
    P1LumpedInvMassOperator invMass;
    bool solveExactly_;
};

}

