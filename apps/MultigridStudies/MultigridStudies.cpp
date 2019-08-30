
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/timing/TimingJSON.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "tinyhhg_core/composites/P2P2StokesFunction.hpp"
#include "tinyhhg_core/composites/P2P2UnstableStokesOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1QuadraticProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScMinResSolver.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/FullMultigridSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/SORSmoother.hpp"
#include "tinyhhg_core/solvers/UzawaSmoother.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesBlockDiagonalPreconditioner.hpp"
#include "tinyhhg_core/solvers/preconditioners/StokesPressureBlockPreconditioner.hpp"

#include "sqlite/SQLite.h"

namespace hyteg {

using walberla::int64_c;
using walberla::math::pi;

#define NEUMANN_PROBLEM 0
#define COLLIDING_FLOW 0
#define CONSTANTA_POISSON 1

#if CONSTANTA_POISSON

std::function< real_t( const hyteg::Point3D& ) > exact;
std::function< real_t( const hyteg::Point3D& ) > rhs;

std::function< real_t( const hyteg::Point3D& ) > exactConstanta2D = []( const hyteg::Point3D& x ) {
   return cos( pi * x[0] ) - sin( 2.0 * pi * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > rhsConstanta2D = []( const hyteg::Point3D& x ) {
   return pi * pi * cos( pi * x[0] ) - 4.0 * pi * pi * sin( 2.0 * pi * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > exactConstanta3D = []( const hyteg::Point3D& x ) {
   return cos( pi * x[0] ) - sin( 2.0 * pi * x[1] ) + cos( 2.0 * pi * x[2] );
};

std::function< real_t( const hyteg::Point3D& ) > rhsConstanta3D = []( const hyteg::Point3D& x ) {
   return pi * pi * cos( pi * x[0] ) - 4.0 * pi * pi * sin( 2.0 * pi * x[1] ) + 4.0 * pi * pi * cos( 2.0 * pi * x[2] );
};

#else // END IF CONSTANTA
#if 1
std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) {
   return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > rhs = []( const hyteg::Point3D& x ) {
   return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
};
#else
std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };

std::function< real_t( const hyteg::Point3D& ) > rhs    = []( const hyteg::Point3D& ) { return 0; };
#endif
#endif // END ELSE CONSTANTA

#if NEUMANN_PROBLEM
std::function< real_t( const hyteg::Point3D& ) > bcU = []( const hyteg::Point3D& x ) {
   if ( std::abs( x[0] + 1 ) < 1e-8 )
   {
      return 1 - x[1] * x[1];
   }
   else
   {
      return 0.0;
   }
};

std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& x ) { return 1 - x[1] * x[1]; };

std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& ) { return 0.0; };
std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& x ) { return -2 * x[0]; };
std::function< real_t( const hyteg::Point3D& ) > rhsU   = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > rhsV   = []( const hyteg::Point3D& x ) { return 0; };

#else

#if COLLIDING_FLOW
std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& x ) { return 20 * x[0] * std::pow( x[1], 3.0 ); };
std::function< real_t( const hyteg::Point3D& ) > bcU    = []( const hyteg::Point3D& x ) { return 20 * x[0] * std::pow( x[1], 3.0 ); };
std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& x ) {
   return 5 * std::pow( x[0], 4.0 ) - 5 * std::pow( x[1], 4.0 );
};
std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& x ) {
   return 60 * std::pow( x[0], 2.0 ) * x[1] - 20 * std::pow( x[1], 3.0 );
};
std::function< real_t( const hyteg::Point3D& ) > rhsU = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > rhsV = []( const hyteg::Point3D& ) { return 0; };

#else
std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& x ) {
   return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > bcU = []( const hyteg::Point3D& x ) {
   return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& x ) {
   return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > bcV = []( const hyteg::Point3D& x ) {
    return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > exactW = []( const hyteg::Point3D& ) {
    return real_c(0);
};
std::function< real_t( const hyteg::Point3D& ) > bcW = []( const hyteg::Point3D& ) {
    return real_c(0);
};
std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& x ) {
   return 2.5 * pi * std::cos( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > rhsU = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > rhsV = []( const hyteg::Point3D& x ) {
   return -12.5 * pi * pi * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > rhsW = []( const hyteg::Point3D& ) {
    return real_c(0);
};
#endif
#endif

std::function< real_t( const hyteg::Point3D& ) > shellExactU = []( const hyteg::Point3D& x ) { return -4 * std::cos( 4 * x[2] ); };
std::function< real_t( const hyteg::Point3D& ) > shellExactV = []( const hyteg::Point3D& x ) { return 8 * std::cos( 8 * x[0] ); };
std::function< real_t( const hyteg::Point3D& ) > shellExactW = []( const hyteg::Point3D& x ) { return -2 * std::cos( 2 * x[1] ); };
std::function< real_t( const hyteg::Point3D& ) > shellExactP = []( const hyteg::Point3D& x ) {
   return std::sin( 4 * x[0] ) * std::sin( 8 * x[1] ) * std::sin( 2 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > shellRhsU = []( const hyteg::Point3D& x ) {
   return 4 * std::sin( 8 * x[1] ) * std::sin( 2 * x[2] ) * std::cos( 4 * x[0] ) - 64 * std::cos( 4 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > shellRhsV = []( const hyteg::Point3D& x ) {
   return 8 * std::sin( 4 * x[0] ) * std::sin( 2 * x[2] ) * std::cos( 8 * x[1] ) + 512 * std::cos( 8 * x[0] );
};
std::function< real_t( const hyteg::Point3D& ) > shellRhsW = []( const hyteg::Point3D& x ) {
   return 2 * std::sin( 4 * x[0] ) * std::sin( 8 * x[1] ) * std::cos( 2 * x[2] ) - 8 * std::cos( 2 * x[1] );
};

template < typename Function, typename LaplaceOperator, typename MassOperator >
void calculateErrorAndResidual( const uint_t&          level,
                                const LaplaceOperator& A,
                                const MassOperator&,
                                const Function& u,
                                const Function& f,
                                const Function& uExact,
                                const Function& error,
                                const Function& residual,
                                const Function& tmp,
                                long double&    LInfError,
                                long double&    L2Error,
                                long double&    LInfResidual,
                                long double&    L2Residual )
{
   error.assign( {1.0, -1.0}, {uExact, u}, level, All );

   tmp.interpolate( real_c( 0 ), level, All );
   A.apply( u, tmp, level, Inner );
   residual.assign( {1.0, -1.0}, {f, tmp}, level, All );

   auto num = numberOfGlobalDoFs< typename Function::Tag >( *u.getStorage(), level );

   LInfError    = error.getMaxMagnitude( level );
   LInfResidual = residual.getMaxMagnitude( level );

   L2Error    = std::sqrt( error.dotGlobal( error, level, Inner ) / (long double) ( num ) );
   L2Residual = std::sqrt( residual.dotGlobal( residual, level, Inner ) / (long double) ( num ) );
}

template < typename Function, typename StokesOperator >
void calculateErrorAndResidualStokes( const uint_t&         level,
                                      const StokesOperator& A,
                                      const Function&       u,
                                      const Function&       f,
                                      const Function&       error,
                                      long double&          l2ErrorU,
                                      long double&          l2ErrorP,
                                      long double&          l2ResidualU,
                                      long double&          l2ResidualP )
{
   auto numU = numberOfGlobalDoFs< typename Function::VelocityFunction_T::Tag >( *u.u.getStorage(), level );
   auto numP = numberOfGlobalDoFs< typename Function::PressureFunction_T::Tag >( *u.p.getStorage(), level );

   // residual (storing in error function to minimize mem overhead)

   error.interpolate( real_c( 0 ), level, All );
   A.apply( u, error, level, Inner | NeumannBoundary );
   error.assign( {1.0, -1.0}, {f, error}, level, All );

   real_t sumVelocityResidualDot = 0.0;
   sumVelocityResidualDot += error.u.dotGlobal( error.u, level, Inner | NeumannBoundary );
   sumVelocityResidualDot += error.v.dotGlobal( error.v, level, Inner | NeumannBoundary );
   if ( error.w.getStorage()->hasGlobalCells() )
   {
      sumVelocityResidualDot += error.w.dotGlobal( error.w, level, Inner | NeumannBoundary );
      sumVelocityResidualDot /= real_c( (long double) ( 3 * numU ) );
   }
   else
   {
      sumVelocityResidualDot /= real_c( (long double) ( 2 * numU ) );
   }

   l2ResidualU = std::sqrt( sumVelocityResidualDot );
   l2ResidualP = std::sqrt( error.p.dotGlobal( error.p, level, Inner | NeumannBoundary ) / (long double) ( numP ) );

   // error

   error.u.interpolate( exactU, level, All );
   error.v.interpolate( exactV, level, All );
   error.w.interpolate( exactW, level, All );
   error.p.interpolate( exactP, level, All );
   error.assign( {1.0, -1.0}, {error, u}, level, All );

   vertexdof::projectMean( error.p, level );

   real_t sumVelocityErrorDot = 0.0;
   sumVelocityErrorDot += error.u.dotGlobal( error.u, level, Inner | NeumannBoundary );
   sumVelocityErrorDot += error.v.dotGlobal( error.v, level, Inner | NeumannBoundary );
   if ( error.w.getStorage()->hasGlobalCells() )
   {
      sumVelocityErrorDot += error.w.dotGlobal( error.w, level, Inner | NeumannBoundary );
      sumVelocityErrorDot /= real_c( (long double) ( 3 * numU ) );
   }
   else
   {
      sumVelocityErrorDot /= real_c( (long double) ( 2 * numU ) );
   }

   l2ErrorU = std::sqrt( sumVelocityErrorDot );
   l2ErrorP = std::sqrt( error.p.dotGlobal( error.p, level, Inner | NeumannBoundary ) / (long double) ( numP ) );
}

template < typename Function, typename LaplaceOperator, typename MassOperator >
void calculateDiscretizationError( const std::shared_ptr< PrimitiveStorage >& storage,
                                   const uint_t&                              level,
                                   long double&                               l2DiscretizationError )
{
   Function u( "u", storage, level, level );
   Function f( "f", storage, level, level );

   Function uExact( "uExact", storage, level, level );
   Function residual( "residual", storage, level, level );
   Function error( "error", storage, level, level );
   Function tmp( "tmp", storage, level, level );

   LaplaceOperator A( storage, level, level );
   MassOperator    M( storage, level, level );

   u.interpolate( exact, level, DirichletBoundary );
   uExact.interpolate( exact, level, All );

   tmp.interpolate( rhs, level, All );
   M.apply( tmp, f, level, All );

#ifdef HHG_BUILD_WITH_PETSC
   auto solver = std::make_shared< PETScLUSolver< LaplaceOperator > >( storage, level );
#else
   auto solver = std::make_shared< CGSolver< LaplaceOperator > >( storage, level, level );
#endif
   solver->solve( A, u, f, level );

   long double L2Error;
   long double l2Residual;
   long double L2Residual;
   calculateErrorAndResidual(
       level, A, M, u, f, uExact, error, residual, tmp, L2Error, l2DiscretizationError, l2Residual, L2Residual );
}

template < typename StokesFunction, typename StokesOperator, typename MassOperator >
void calculateDiscretizationErrorStokes( const std::shared_ptr< PrimitiveStorage >& storage,
                                         const uint_t&                              level,
                                         long double&                               l2DiscretizationErrorU,
                                         long double&                               l2DiscretizationErrorP )
{
   StokesFunction u( "u", storage, level, level );
   StokesFunction f( "f", storage, level, level );

   StokesFunction error( "error", storage, level, level );
   StokesFunction tmp( "tmp", storage, level, level );

   StokesOperator A( storage, level, level );
   MassOperator   M( storage, level, level );

   u.u.interpolate( bcU, level, DirichletBoundary );
   u.v.interpolate( bcV, level, DirichletBoundary );
   u.w.interpolate( bcW, level, DirichletBoundary );

   tmp.u.interpolate( rhsU, level, All );
   tmp.v.interpolate( rhsV, level, All );
   tmp.w.interpolate( rhsW, level, All );
   M.apply( tmp.u, f.u, level, All );
   M.apply( tmp.v, f.v, level, All );
   M.apply( tmp.w, f.w, level, All );

#ifdef HHG_BUILD_WITH_PETSC
   auto solver = std::make_shared< PETScMinResSolver< StokesOperator > >( storage, level, 1e-16 );
   // auto solver = std::make_shared< PETScLUSolver< StokesOperator > >( storage, level );
#else
   auto cgVelocity =
       std::make_shared< CGSolver< typename StokesOperator::VelocityOperator_T > >( storage, level, level, 0, 1e-14 );
   auto preconditioner = std::make_shared< StokesPressureBlockPreconditioner< StokesOperator, P1LumpedInvMassOperator > >(
       storage, level, level ); //, 1, cgVelocity );
   auto solver = std::make_shared< MinResSolver< StokesOperator > >(
       storage, level, level, std::numeric_limits< uint_t >::max(), 1e-16, preconditioner );
#endif
   solver->solve( A, u, f, level );

   vertexdof::projectMean( u.p, level );

   long double l2ResidualU;
   long double l2ResidualP;

   calculateErrorAndResidualStokes(
       level, A, u, f, error, l2DiscretizationErrorU, l2DiscretizationErrorP, l2ResidualU, l2ResidualP );
}

template < typename Function,
           typename LaplaceOperator,
           typename MassOperator,
           typename Restriction,
           typename Prolongation,
           typename FMGProlongation >
void MultigridLaplace( const std::shared_ptr< PrimitiveStorage >&           storage,
                       const uint_t&                                        minLevel,
                       const uint_t&                                        maxLevel,
                       const uint_t&                                        numCycles,
                       const CycleType                                      cycleType,
                       const uint_t&                                        fmgInnerCycles,
                       const real_t&                                        L2residualTolerance,
                       const real_t&                                        sorRelax,
                       const uint_t&                                        preSmoothingSteps,
                       const uint_t&                                        postSmoothingSteps,
                       const bool&                                          outputVTK,
                       const uint_t&                                        skipCyclesForAvgConvRate,
                       const bool&                                          calcDiscretizationError,
                       std::map< std::string, walberla::int64_t >&          sqlIntegerProperties,
                       std::map< std::string, double >&                     sqlRealProperties,
                       std::map< std::string, std::string >&                sqlStringProperties,
                       std::map< uint_t, std::map< std::string, double > >& sqlRealPropertiesMG )
{
   WALBERLA_UNUSED( sqlStringProperties );

   Function u( "u", storage, minLevel, maxLevel );
   Function f( "f", storage, minLevel, maxLevel );

   Function uExact( "uExact", storage, minLevel, maxLevel );
   Function residual( "residual", storage, minLevel, maxLevel );
   Function error( "error", storage, minLevel, maxLevel );
   Function tmp( "tmp", storage, minLevel, maxLevel );

   LaplaceOperator A( storage, minLevel, maxLevel );
   MassOperator    M( storage, minLevel, maxLevel );

   long double LInfError;
   long double L2Error;
   long double LInfResidual;
   long double L2Residual;

   ////////////////////
   // Initialize VTK //
   ////////////////////

   VTKOutput vtkOutput( "vtk", "P2MultigridLaplace", storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( uExact );
   vtkOutput.add( residual );
   vtkOutput.add( error );

   //////////////////////////////////////////////
   // Initialize functions and right-hand side //
   //////////////////////////////////////////////

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      u.interpolate( exact, level, DirichletBoundary );
      uExact.interpolate( exact, level, All );

      tmp.interpolate( rhs, level, All );
      M.apply( tmp, f, level, All );
   }

   /////////////////////////
   // Misc setup and info //
   /////////////////////////

   long double discretizationError = 0.0;
   if ( calcDiscretizationError )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "l2 discretization error per level:" );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         calculateDiscretizationError< Function, LaplaceOperator, MassOperator >( storage, level, discretizationError );
         WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::scientific << discretizationError );
         sqlRealProperties["l2_discr_error_level_" + std::to_string( level )] = real_c( discretizationError );
      }
      WALBERLA_LOG_DEVEL_ON_ROOT( "" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Number of unknowns (including boundary):" )
   uint_t totalDoFs = 0;
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      const uint_t dofsThisLevel = numberOfGlobalDoFs< typename Function::Tag >( *storage, level );
      WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::setw( 15 ) << dofsThisLevel );
      totalDoFs += dofsThisLevel;
   }
   WALBERLA_LOG_INFO_ON_ROOT( " ----------------------------- " );
   WALBERLA_LOG_INFO_ON_ROOT( "  total:    " << std::setw( 15 ) << totalDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   sqlIntegerProperties["total_dofs"] = int64_c( totalDoFs );

   walberla::WcTimer timer;
   double            timeError;
   double            timeVTK;
   double            timeCycle;

   timer.reset();
   calculateErrorAndResidual( maxLevel, A, M, u, f, uExact, error, residual, tmp, LInfError, L2Error, LInfResidual, L2Residual );
   timer.end();
   timeError = timer.last();

   timer.reset();
   if ( outputVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }
   timer.end();
   timeVTK = timer.last();

   WALBERLA_LOG_INFO_ON_ROOT(
       " After cycle... ||   LInf error |     L2 error | L2 error reduction ||  LInf res.   |  L2 residual | L2 residual reduction || time cycle [s] | time error calculation [s] | time VTK [s] |" );
   WALBERLA_LOG_INFO_ON_ROOT(
       " ---------------++--------------+--------------+--------------------++--------------+--------------+-----------------------++----------------+----------------------------+--------------|" );
   WALBERLA_LOG_INFO_ON_ROOT( "        initial || " << std::scientific << LInfError << " | " << L2Error << " | "
                                                    << "               --- || " << LInfResidual << " | " << L2Residual
                                                    << " |                   --- ||            --- | " << std::fixed
                                                    << std::setprecision( 2 ) << std::setw( 26 ) << timeError << " | "
                                                    << std::setw( 12 ) << timeVTK << " |" );

   long double avgL2ErrorConvergenceRate      = 0;
   long double avgL2ResidualConvergenceRate   = 0;
   long double avgLInfErrorConvergenceRate    = 0;
   long double avgLInfResidualConvergenceRate = 0;

   long double L2ErrorReduction      = 0;
   long double L2ResidualReduction   = 0;
   long double LInfErrorReduction    = 0;
   long double LInfResidualReduction = 0;

   sqlRealPropertiesMG[0]["L2_error"]              = real_c( L2Error );
   sqlRealPropertiesMG[0]["L2_error_reduction"]    = real_c( L2ErrorReduction );
   sqlRealPropertiesMG[0]["L2_residual"]           = real_c( L2Residual );
   sqlRealPropertiesMG[0]["L2_residual_reduction"] = real_c( L2ResidualReduction );

   sqlRealPropertiesMG[0]["LInf_error"]              = real_c( LInfError );
   sqlRealPropertiesMG[0]["LInf_error_reduction"]    = real_c( LInfErrorReduction );
   sqlRealPropertiesMG[0]["LInf_residual"]           = real_c( LInfResidual );
   sqlRealPropertiesMG[0]["LInf_residual_reduction"] = real_c( LInfResidualReduction );

   ///////////
   // Solve //
   ///////////

   auto smoother = std::make_shared< SORSmoother< LaplaceOperator > >( sorRelax );
#ifdef HHG_BUILD_WITH_PETSC
   auto coarseGridSolver = std::make_shared< PETScLUSolver< LaplaceOperator > >( storage, minLevel );
#else
   auto coarseGridSolver = std::make_shared< CGSolver< LaplaceOperator > >( storage, minLevel, minLevel );
#endif

   auto prolongationOperator = std::make_shared< Prolongation >();
   auto restrictionOperator  = std::make_shared< Restriction >();

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< LaplaceOperator > >( storage,
                                                                                           smoother,
                                                                                           coarseGridSolver,
                                                                                           restrictionOperator,
                                                                                           prolongationOperator,
                                                                                           minLevel,
                                                                                           maxLevel,
                                                                                           preSmoothingSteps,
                                                                                           postSmoothingSteps,
                                                                                           0,
                                                                                           cycleType );

   auto fmgProlongation = std::make_shared< FMGProlongation >();

   auto postCycle = [&]( uint_t currentLevel ) {
      long double _l2Error, _L2Error, _l2Residual, _L2Residual;
      calculateErrorAndResidual(
          currentLevel, A, M, u, f, uExact, error, residual, tmp, _l2Error, _L2Error, _l2Residual, _L2Residual );
      sqlRealProperties["fmg_l2_error_level_" + std::to_string( currentLevel )] = real_c( _l2Error );
      WALBERLA_LOG_INFO_ON_ROOT( "    fmg level " << currentLevel << ": l2 error: " << std::scientific << _l2Error );
   };

   FullMultigridSolver< LaplaceOperator > fullMultigridSolver(
       storage, multigridSolver, fmgProlongation, minLevel, maxLevel, fmgInnerCycles, postCycle );

   uint_t numExecutedCycles = 0;
   for ( uint_t cycle = 1; cycle <= numCycles; cycle++ )
   {
      const long double lastL2Error    = L2Error;
      const long double lastL2Residual = L2Residual;

      const long double lastLInfError    = LInfError;
      const long double lastLInfResidual = LInfResidual;

      timer.reset();
      if ( cycle == 1 && fmgInnerCycles > 0 )
      {
         fullMultigridSolver.solve( A, u, f, maxLevel );
      }
      else
      {
         multigridSolver->solve( A, u, f, maxLevel );
      }
      timer.end();
      timeCycle = timer.last();

      numExecutedCycles++;

      timer.reset();
      calculateErrorAndResidual(
          maxLevel, A, M, u, f, uExact, error, residual, tmp, LInfError, L2Error, LInfResidual, L2Residual );
      timer.end();
      timeError = timer.last();

      timer.reset();
      if ( outputVTK )
      {
         vtkOutput.write( maxLevel, cycle );
      }
      timer.end();
      timeVTK = timer.last();

      L2ErrorReduction      = L2Error / lastL2Error;
      L2ResidualReduction   = L2Residual / lastL2Residual;
      LInfErrorReduction    = LInfError / lastLInfError;
      LInfResidualReduction = LInfResidual / lastLInfResidual;

      WALBERLA_LOG_INFO_ON_ROOT(
          std::setw( 15 ) << cycle << " || " << std::scientific << LInfError << " | " << L2Error << " | "
                          << "      " << L2ErrorReduction << " || " << LInfResidual << " | " << L2Residual << " |          "
                          << L2ResidualReduction << " || " << std::fixed << std::setprecision( 2 ) << std::setw( 14 ) << timeCycle
                          << " | " << std::setw( 26 ) << timeError << " | " << std::setw( 12 ) << timeVTK << " | "
                          << " | ratio discr.err: " << ( calcDiscretizationError ? LInfError / discretizationError : 0.0 ) );

      if ( cycle > skipCyclesForAvgConvRate )
      {
         avgL2ErrorConvergenceRate += L2ErrorReduction;
         avgL2ResidualConvergenceRate += L2ResidualReduction;
         avgLInfErrorConvergenceRate += LInfErrorReduction;
         avgLInfResidualConvergenceRate += LInfResidualReduction;
      }

      sqlRealPropertiesMG[cycle]["L2_error"]              = real_c( L2Error );
      sqlRealPropertiesMG[cycle]["L2_error_reduction"]    = real_c( L2ErrorReduction );
      sqlRealPropertiesMG[cycle]["L2_residual"]           = real_c( L2Residual );
      sqlRealPropertiesMG[cycle]["L2_residual_reduction"] = real_c( L2ResidualReduction );

      sqlRealPropertiesMG[cycle]["LInf_error"]              = real_c( LInfError );
      sqlRealPropertiesMG[cycle]["LInf_error_reduction"]    = real_c( LInfErrorReduction );
      sqlRealPropertiesMG[cycle]["LInf_residual"]           = real_c( LInfResidual );
      sqlRealPropertiesMG[cycle]["LInf_residual_reduction"] = real_c( LInfResidualReduction );

      if ( L2Residual < L2residualTolerance )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "L2 residual dropped below tolerance." )
         break;
      }
   }

   avgL2ErrorConvergenceRate /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );
   avgL2ResidualConvergenceRate /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );

   avgLInfErrorConvergenceRate /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );
   avgLInfResidualConvergenceRate /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );

   sqlRealProperties["avg_capital_L2_error_conv_rate"]    = real_c( avgL2ErrorConvergenceRate );
   sqlRealProperties["avg_capital_L2_residual_conv_rate"] = real_c( avgL2ResidualConvergenceRate );

   sqlRealProperties["avg_LInf_error_conv_rate"]    = real_c( avgLInfErrorConvergenceRate );
   sqlRealProperties["avg_LInf_residual_conv_rate"] = real_c( avgLInfResidualConvergenceRate );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Average convergence rates:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 error:      " << std::scientific << avgL2ErrorConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual:   " << std::scientific << avgL2ResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - LInf error:    " << std::scientific << avgLInfErrorConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - LInf residual: " << std::scientific << avgLInfResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

template < typename StokesOperator, typename StokesFunction >
void DCStokesRHSSetup( const std::shared_ptr< PrimitiveStorage >&,
                       const uint_t&,
                       const StokesOperator&,
                       const StokesFunction&,
                       const P1StokesFunction< real_t >& )
{
   WALBERLA_ABORT( "Defect correction not implemented for this discretization." )
}

template <>
void DCStokesRHSSetup< P1StokesOperator, P1StokesFunction< real_t > >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                       const uint_t&                              p1Level,
                                                                       const P1StokesOperator&           p1StokesOperator,
                                                                       const P1StokesFunction< real_t >& u,
                                                                       const P1StokesFunction< real_t >& p1DefectCorrectionRHS )
{
   const uint_t p2Level = p1Level - 1;

   walberla::WcTimer timer;
   timer.reset();

   P1StokesFunction< real_t > Au_P1( "Au_P1", storage, p1Level, p1Level );
   P1StokesFunction< real_t > Au_P2_converted_to_P1( "Au_P1", storage, p1Level, p1Level );
   P1StokesFunction< real_t > f_P2_on_P1_space( "f_P2_on_p1_space", storage, p1Level, p1Level );

   P2P2StokesFunction< real_t > u_P2( "u_P2", storage, p2Level, p2Level );
   P2P2StokesFunction< real_t > Au_P2( "Au_P2", storage, p2Level, p2Level );
   P2P2StokesFunction< real_t > tmp_P2( "tmp_P2", storage, p2Level, p2Level );
   P2P2StokesFunction< real_t > f_P2( "f_P2", storage, p2Level, p2Level );

   P2P2UnstableStokesOperator A_P2( storage, p2Level, p2Level );
   P2ConstantMassOperator     M_P2( storage, p2Level, p2Level );

   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "-> time DC: function allocation: " << std::fixed << std::setprecision( 2 ) << std::setw( 10 )
                                                                  << timer.last() << "s" )

   timer.reset();
   // set up higher order RHS
   tmp_P2.u.interpolate( rhsU, p2Level, All );
   tmp_P2.v.interpolate( rhsV, p2Level, All );
   M_P2.apply( tmp_P2.u, f_P2.u, p2Level, All );
   M_P2.apply( tmp_P2.v, f_P2.v, p2Level, All );
   f_P2_on_P1_space.u.assign( f_P2.u, p1Level, All );
   f_P2_on_P1_space.v.assign( f_P2.v, p1Level, All );

   // A * u (linear)
   p1StokesOperator.apply( u, Au_P1, p1Level, Inner );

   // A_higher_order * u (quadratic)
   // u_quadratic is given by direct injection of the linear coefficients
   u_P2.u.assign( u.u, p2Level, All );
   u_P2.v.assign( u.v, p2Level, All );
   u_P2.w.assign( u.w, p2Level, All );
   u_P2.p.assign( u.p, p2Level, All );

   A_P2.apply( u_P2, Au_P2, p2Level, Inner );

   Au_P2_converted_to_P1.u.assign( Au_P2.u, p1Level, All );
   Au_P2_converted_to_P1.v.assign( Au_P2.v, p1Level, All );
   Au_P2_converted_to_P1.w.assign( Au_P2.w, p1Level, All );
   Au_P2_converted_to_P1.p.assign( Au_P2.p, p1Level, All );

   // defect correction
   // f_correction = f - (A_higher_order * u^i-1) + (A * u^i-1)
   p1DefectCorrectionRHS.assign( {1.0, -1.0, 1.0}, {f_P2_on_P1_space, Au_P2_converted_to_P1, Au_P1}, p1Level, All );
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "-> time DC: RHS calculation:     " << std::fixed << std::setprecision( 2 ) << std::setw( 10 )
                                                                  << timer.last() << "s" )
}

template < typename StokesOperator, typename StokesFunction >
void DCStokesRunCycle( const std::shared_ptr< GeometricMultigridSolver< StokesOperator > >&,
                       const StokesOperator&,
                       const StokesFunction&,
                       const P1StokesFunction< real_t >&,
                       const uint_t& )
{
   WALBERLA_ABORT( "Defect correction not implemented for this discretization." )
}

template <>
void DCStokesRunCycle< P1StokesOperator, P1StokesFunction< real_t > >(
    const std::shared_ptr< GeometricMultigridSolver< P1StokesOperator > >& solver,
    const P1StokesOperator&                                                p1StokesOperator,
    const P1StokesFunction< real_t >&                                      u,
    const P1StokesFunction< real_t >&                                      f_dc,
    const uint_t&                                                          level )
{
   solver->solve( p1StokesOperator, u, f_dc, level );
}

template < typename StokesFunction,
           typename StokesFunctionNumerator,
           typename StokesOperator,
           typename MassOperator,
           typename Restriction,
           typename Prolongation,
           typename FMGProlongation >
void MultigridStokes( const std::shared_ptr< PrimitiveStorage >&           storage,
                      const uint_t&                                        minLevel,
                      const uint_t&                                        maxLevel,
                      const uint_t&                                        numCycles,
                      const CycleType                                      cycleType,
                      const uint_t&                                        fmgInnerCycles,
                      const real_t&                                        L2residualTolerance,
                      const real_t&                                        sorRelax,
                      const bool&                                          symmGSVelocity,
                      const uint_t&                                        numGSVelocity,
                      const bool&                                          symmGSPressure,
                      const uint_t&                                        numGSPressure,
                      const uint_t&                                        preSmoothingSteps,
                      const uint_t&                                        postSmoothingSteps,
                      const uint_t&                                        smoothingIncrement,
                      const bool&                                          projectPressureAfterRestriction,
                      const uint_t&                                        coarseGridMaxIterations,
                      const real_t&                                        coarseResidualTolerance,
                      const bool&                                          outputVTK,
                      const uint_t&                                        skipCyclesForAvgConvRate,
                      const bool&                                          calcDiscretizationError,
                      const uint_t&                                        cyclesBeforeDC,
                      const uint_t&                                        postDCPreSmoothingSteps,
                      const uint_t&                                        postDCPostSmoothingSteps,
                      const uint_t&                                        postDCSmoothingIncrement,
                      std::map< std::string, walberla::int64_t >&          sqlIntegerProperties,
                      std::map< std::string, double >&                     sqlRealProperties,
                      std::map< std::string, std::string >&                sqlStringProperties,
                      std::map< uint_t, std::map< std::string, double > >& sqlRealPropertiesMG )
{
   WALBERLA_UNUSED( sqlStringProperties );

   if ( cyclesBeforeDC > 0 )
   {
      if ( !std::is_same< typename StokesFunction::Tag, P1StokesFunctionTag >::value )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "DC enabled, but only works with P1-P1-stab discretization!" )
      }
   }

   StokesFunction u( "u", storage, minLevel, maxLevel );
   StokesFunction f( "f", storage, minLevel, maxLevel );

   StokesFunction error( "error", storage, minLevel, maxLevel );

   std::shared_ptr< P1StokesFunction< real_t > > f_dc;
   if ( cyclesBeforeDC > 0 )
      f_dc = std::make_shared< P1StokesFunction< real_t > >( "f_dc", storage, minLevel, maxLevel );

   StokesOperator A( storage, minLevel, maxLevel );
   MassOperator   M( storage, minLevel, maxLevel );

   long double l2ErrorU;
   long double l2ErrorP;
   long double l2ResidualU;
   long double l2ResidualP;

   ////////////////////
   // Initialize VTK //
   ////////////////////

   VTKOutput vtkOutput( "vtk", "MultigridStudies", storage );
   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( error );

   ///////////////////////////
   // Output exact solution //
   ///////////////////////////

   if ( outputVTK )
   {
      error.u.interpolate( exactU, maxLevel );
      error.v.interpolate( exactV, maxLevel );
      error.w.interpolate( exactW, maxLevel );
      error.p.interpolate( exactP, maxLevel );

      VTKOutput exactSolutionVTKOutput( "vtk", "MultigridStudiesExact", storage );
      exactSolutionVTKOutput.add( error );
      exactSolutionVTKOutput.write( maxLevel );
   }

   //////////////////////////////////////////////
   // Initialize functions and right-hand side //
   //////////////////////////////////////////////

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      u.u.interpolate( bcU, level, DirichletBoundary );
      u.v.interpolate( bcV, level, DirichletBoundary );
      u.w.interpolate( bcW, level, DirichletBoundary );

      // using error as tmp function here
      error.u.interpolate( rhsU, level, All );
      error.v.interpolate( rhsV, level, All );
      error.w.interpolate( rhsW, level, All );
      M.apply( error.u, f.u, level, All );
      M.apply( error.v, f.v, level, All );
      M.apply( error.w, f.w, level, All );
   }

   /////////////////////////
   // Misc setup and info //
   /////////////////////////

   long double discretizationErrorU = 0;
   long double discretizationErrorP = 0;
   if ( calcDiscretizationError )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "l2 discretization error ( u | p ) per level:" );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         calculateDiscretizationErrorStokes< StokesFunction, StokesOperator, MassOperator >(
             storage, level, discretizationErrorU, discretizationErrorP );
         WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::scientific << discretizationErrorU
                                               << " | " << discretizationErrorP );
         sqlRealProperties["l2_discr_error_u_level_" + std::to_string( level )] = real_c( discretizationErrorU );
         sqlRealProperties["l2_discr_error_p_level_" + std::to_string( level )] = real_c( discretizationErrorP );
      }
      WALBERLA_LOG_DEVEL_ON_ROOT( "" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Number of unknowns (excluding boundary for velocity):" )
   uint_t totalDoFs = 0;
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      const uint_t velocityDoFsThisLevel =
          storage->hasGlobalCells() ?
              3 * numberOfGlobalInnerDoFs< typename StokesFunction::VelocityFunction_T::Tag >( *storage, level ) :
              2 * numberOfGlobalInnerDoFs< typename StokesFunction::VelocityFunction_T::Tag >( *storage, level );
      const uint_t dofsThisLevel =
          numberOfGlobalDoFs< typename StokesFunction::PressureFunction_T::Tag >( *storage, level ) + velocityDoFsThisLevel;
      WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::setw( 15 ) << dofsThisLevel );
      sqlIntegerProperties["total_dofs_level_" + std::to_string( level )] = int64_c( dofsThisLevel );
      totalDoFs += dofsThisLevel;
   }
   WALBERLA_LOG_INFO_ON_ROOT( " ----------------------------- " );
   WALBERLA_LOG_INFO_ON_ROOT( "  total:    " << std::setw( 15 ) << totalDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   sqlIntegerProperties["total_dofs"] = int64_c( totalDoFs );

   storage->getTimingTree()->reset();

   walberla::WcTimer timer;
   walberla::WcTimer timerFMGErrorCalculation;
   double            timeError;
   double            timeVTK;
   double            timeCycle;

   ///////////
   // Solve //
   ///////////

   const uint_t coarseGridMaxLevel = ( numCycles == 0 ? maxLevel : minLevel );

   auto smoother = std::make_shared< UzawaSmoother< StokesOperator > >( storage,
                                                                        error,
                                                                        minLevel,
                                                                        maxLevel,
                                                                        sorRelax,
                                                                        Inner | NeumannBoundary,
                                                                        symmGSVelocity,
                                                                        numGSVelocity,
                                                                        symmGSPressure,
                                                                        numGSPressure );

#ifdef HHG_BUILD_WITH_PETSC
   //   auto petscSolver = std::make_shared< PETScMinResSolver< StokesOperator > >(
   //       storage, coarseGridMaxLevel, coarseResidualTolerance, coarseGridMaxIterations );
   auto petscSolver = std::make_shared< PETScLUSolver< StokesOperator > >( storage, coarseGridMaxLevel );
   WALBERLA_UNUSED( coarseGridMaxIterations );
   WALBERLA_UNUSED( coarseResidualTolerance );
#else
   //   const uint_t preconditionerCGIterations = 0;
   //
   //   auto cgVelocity = std::make_shared< CGSolver< typename StokesOperator::VelocityOperator_T > >(
   //       storage, minLevel, coarseGridMaxLevel, preconditionerCGIterations, 1e-14 );

   auto preconditioner = std::make_shared< StokesPressureBlockPreconditioner< StokesOperator, P1LumpedInvMassOperator > >(
       storage, minLevel, coarseGridMaxLevel ); //, 1, cgVelocity );
   auto coarseGridSolver = std::make_shared< MinResSolver< StokesOperator > >(
       storage, minLevel, coarseGridMaxLevel, coarseGridMaxIterations, coarseResidualTolerance, preconditioner );
#endif

   auto prolongationOperator = std::make_shared< Prolongation >();
   auto restrictionOperator  = std::make_shared< Restriction >( projectPressureAfterRestriction );

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >( storage,
                                                                                          error,
                                                                                          smoother,
#ifdef HHG_BUILD_WITH_PETSC
                                                                                          petscSolver,
#else
                                                                                          coarseGridSolver,
#endif
                                                                                          restrictionOperator,
                                                                                          prolongationOperator,
                                                                                          minLevel,
                                                                                          maxLevel,
                                                                                          preSmoothingSteps,
                                                                                          postSmoothingSteps,
                                                                                          smoothingIncrement,
                                                                                          cycleType );

   auto fmgProlongation = std::make_shared< FMGProlongation >();

   auto postCycle = [&]( uint_t currentLevel ) {
      timerFMGErrorCalculation.start();
      long double _l2ErrorU, _l2ErrorP, _l2ResidualU, _l2ResidualP;
      calculateErrorAndResidualStokes( currentLevel, A, u, f, error, _l2ErrorU, _l2ErrorP, _l2ResidualU, _l2ResidualP );
      sqlRealProperties["fmg_l2_error_u_level_" + std::to_string( currentLevel )]    = real_c( _l2ErrorU );
      sqlRealProperties["fmg_l2_error_p_level_" + std::to_string( currentLevel )]    = real_c( _l2ErrorP );
      sqlRealProperties["fmg_l2_residual_u_level_" + std::to_string( currentLevel )] = real_c( _l2ErrorU );
      sqlRealProperties["fmg_l2_residual_p_level_" + std::to_string( currentLevel )] = real_c( _l2ErrorP );
      timerFMGErrorCalculation.end();
      WALBERLA_LOG_INFO_ON_ROOT( "    fmg level " << currentLevel << ": l2 error u: " << std::scientific << _l2ErrorU
                                                  << " / l2 error p: " << std::scientific << _l2ErrorP << std::fixed
                                                  << " - time error calc: " << timerFMGErrorCalculation.last() << " sec" );
   };

   FullMultigridSolver< StokesOperator > fullMultigridSolver(
       storage, multigridSolver, fmgProlongation, minLevel, maxLevel, fmgInnerCycles, postCycle );

   printFunctionAllocationInfo( *storage, 1 );

   timer.reset();
   calculateErrorAndResidualStokes( maxLevel, A, u, f, error, l2ErrorU, l2ErrorP, l2ResidualU, l2ResidualP );
   timer.end();
   timeError = timer.last();

   timer.reset();
   if ( outputVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }
   timer.end();
   timeVTK = timer.last();

   WALBERLA_LOG_INFO_ON_ROOT(
       " After cycle... ||   l2 error u |   l2 error p |     l2 error u red || l2 residualU | l2 residualP |     l2 residual u red || time cycle [s] | time error calculation [s] | time VTK [s] |" );
   WALBERLA_LOG_INFO_ON_ROOT(
       " ---------------++--------------+--------------+--------------------++--------------+--------------+-----------------------++----------------+----------------------------+--------------|" );
   WALBERLA_LOG_INFO_ON_ROOT( "        initial || " << std::scientific << l2ErrorU << " | " << l2ErrorP << " | "
                                                    << "               --- || " << l2ResidualU << " | " << l2ResidualP
                                                    << " |                   --- ||            --- | " << std::fixed
                                                    << std::setprecision( 2 ) << std::setw( 26 ) << timeError << " | "
                                                    << std::setw( 12 ) << timeVTK << " |" );

   long double avgl2ErrorConvergenceRateU    = 0;
   long double avgl2ResidualConvergenceRateU = 0;
   long double avgl2ErrorConvergenceRateP    = 0;
   long double avgl2ResidualConvergenceRateP = 0;

   long double l2ErrorReductionU    = 0;
   long double l2ResidualReductionU = 0;
   long double l2ErrorReductionP    = 0;
   long double l2ResidualReductionP = 0;

   sqlRealPropertiesMG[0]["l2_error_u"]              = real_c( l2ErrorU );
   sqlRealPropertiesMG[0]["l2_error_reduction_u"]    = real_c( l2ErrorReductionU );
   sqlRealPropertiesMG[0]["l2_residual_u"]           = real_c( l2ResidualU );
   sqlRealPropertiesMG[0]["l2_residual_reduction_u"] = real_c( l2ResidualReductionU );

   sqlRealPropertiesMG[0]["l2_error_p"]              = real_c( l2ErrorP );
   sqlRealPropertiesMG[0]["l2_error_reduction_p"]    = real_c( l2ErrorReductionP );
   sqlRealPropertiesMG[0]["l2_residual_p"]           = real_c( l2ResidualP );
   sqlRealPropertiesMG[0]["l2_residual_reduction_p"] = real_c( l2ResidualReductionP );

   if ( numCycles == 0 )
   {
      timer.reset();
#ifdef HHG_BUILD_WITH_PETSC
      petscSolver->solve( A, u, f, maxLevel );
#else
      coarseGridSolver->solve( A, u, f, maxLevel );
#endif
      timer.end();
      timeCycle = timer.last();
      vertexdof::projectMean( u.p, maxLevel );
      calculateErrorAndResidualStokes( maxLevel, A, u, f, error, l2ErrorU, l2ErrorP, l2ResidualU, l2ResidualP );
      vtkOutput.write( maxLevel, 1 );
      WALBERLA_LOG_INFO_ON_ROOT( std::setw( 15 ) << 1 << " || " << std::scientific << l2ErrorU << " | " << l2ErrorP << " | "
                                                 << "      " << l2ErrorReductionU << " || " << l2ResidualU << " | " << l2ResidualP
                                                 << " |          " << l2ResidualReductionU << " || " << std::fixed
                                                 << std::setprecision( 2 ) << std::setw( 14 ) << timeCycle << " | "
                                                 << std::setw( 26 ) << timeError << " | " << std::setw( 12 ) << timeVTK << " |" );
   }

   uint_t numExecutedCycles = 0;
   for ( uint_t cycle = 1; cycle <= numCycles; cycle++ )
   {
      const long double lastl2ErrorU    = l2ErrorU;
      const long double lastl2ResidualU = l2ResidualU;

      const long double lastl2ErrorP    = l2ErrorP;
      const long double lastl2ResidualP = l2ResidualP;

      if ( cyclesBeforeDC > 0 && numExecutedCycles == cyclesBeforeDC )
      {
         // set up DC RHS once right after the exact number of cycles are performed on the original RHS
         WALBERLA_LOG_INFO_ON_ROOT( "Preparing RHS for DC..." )
         timer.reset();
         DCStokesRHSSetup( storage, maxLevel, A, u, *f_dc );
         timer.end();
         auto timeDCSetup                   = timer.last();
         sqlRealProperties["dc_setup_time"] = timeDCSetup;
         multigridSolver->setSmoothingSteps( postDCPreSmoothingSteps, postDCPostSmoothingSteps, postDCSmoothingIncrement );
      }

      timer.reset();
      if ( cycle == 1 && fmgInnerCycles > 0 )
      {
         LIKWID_MARKER_START( "FMG" );
         fullMultigridSolver.solve( A, u, f, maxLevel );
         LIKWID_MARKER_STOP( "FMG" );
      }
      else if ( cyclesBeforeDC > 0 && numExecutedCycles >= cyclesBeforeDC )
      {
         DCStokesRunCycle( multigridSolver, A, u, *f_dc, maxLevel );
      }
      else
      {
         LIKWID_MARKER_START( "VCYCLE" );
         multigridSolver->solve( A, u, f, maxLevel );
         LIKWID_MARKER_STOP( "VCYCLE" );
      }
      timer.end();
      timeCycle = timer.last();
      if ( cycle == 1 && fmgInnerCycles > 0 )
      {
         timeCycle -= timerFMGErrorCalculation.total();
      }

      numExecutedCycles++;

      if ( !NEUMANN_PROBLEM )
         vertexdof::projectMean( u.p, maxLevel );

      timer.reset();
      calculateErrorAndResidualStokes( maxLevel, A, u, f, error, l2ErrorU, l2ErrorP, l2ResidualU, l2ResidualP );
      timer.end();
      timeError = timer.last();
      if ( cycle == 1 && fmgInnerCycles > 0 )
      {
         timeError += timerFMGErrorCalculation.total();
      }

      timer.reset();
      if ( outputVTK )
      {
         vtkOutput.write( maxLevel, cycle );
      }
      timer.end();
      timeVTK = timer.last();

      l2ErrorReductionU    = l2ErrorU / lastl2ErrorU;
      l2ResidualReductionU = l2ResidualU / lastl2ResidualU;
      l2ErrorReductionP    = l2ErrorP / lastl2ErrorP;
      l2ResidualReductionP = l2ResidualP / lastl2ResidualP;

      WALBERLA_LOG_INFO_ON_ROOT(
          std::setw( 15 ) << cycle << " || " << std::scientific << l2ErrorU << " | " << l2ErrorP << " | "
                          << "      " << l2ErrorReductionU << " || " << l2ResidualU << " | " << l2ResidualP << " |          "
                          << l2ResidualReductionU << " || " << std::fixed << std::setprecision( 2 ) << std::setw( 14 )
                          << timeCycle << " | " << std::setw( 26 ) << timeError << " | " << std::setw( 12 ) << timeVTK
                          << " | ratio discr.err: " << ( calcDiscretizationError ? l2ErrorU / discretizationErrorU : 0.0 ) );

      if ( cycle > skipCyclesForAvgConvRate )
      {
         avgl2ErrorConvergenceRateU += l2ErrorReductionU;
         avgl2ResidualConvergenceRateU += l2ResidualReductionU;
         avgl2ErrorConvergenceRateP += l2ErrorReductionP;
         avgl2ResidualConvergenceRateP += l2ResidualReductionP;
      }

      sqlRealPropertiesMG[cycle]["l2_error_u"]              = real_c( l2ErrorU );
      sqlRealPropertiesMG[cycle]["l2_error_reduction_u"]    = real_c( l2ErrorReductionU );
      sqlRealPropertiesMG[cycle]["l2_residual_u"]           = real_c( l2ResidualU );
      sqlRealPropertiesMG[cycle]["l2_residual_reduction_u"] = real_c( l2ResidualReductionU );

      sqlRealPropertiesMG[cycle]["l2_error_p"]              = real_c( l2ErrorP );
      sqlRealPropertiesMG[cycle]["l2_error_reduction_p"]    = real_c( l2ErrorReductionP );
      sqlRealPropertiesMG[cycle]["l2_residual_p"]           = real_c( l2ResidualP );
      sqlRealPropertiesMG[cycle]["l2_residual_reduction_p"] = real_c( l2ResidualReductionP );

      if ( l2ResidualU < L2residualTolerance )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "l2 residual (u) dropped below tolerance." )
         break;
      }
   }

   avgl2ErrorConvergenceRateU /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );
   avgl2ResidualConvergenceRateU /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );

   avgl2ErrorConvergenceRateP /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );
   avgl2ResidualConvergenceRateP /= real_c( numExecutedCycles - skipCyclesForAvgConvRate );

   sqlRealProperties["avg_l2_error_conv_rate_u"]    = real_c( avgl2ErrorConvergenceRateU );
   sqlRealProperties["avg_l2_residual_conv_rate_u"] = real_c( avgl2ResidualConvergenceRateU );

   sqlRealProperties["avg_l2_error_conv_rate_p"]    = real_c( avgl2ErrorConvergenceRateP );
   sqlRealProperties["avg_l2_residual_conv_rate_p"] = real_c( avgl2ResidualConvergenceRateP );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Average convergence rates:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - l2 error u:    " << std::scientific << avgl2ErrorConvergenceRateU );
   WALBERLA_LOG_INFO_ON_ROOT( "  - l2 residual u: " << std::scientific << avgl2ResidualConvergenceRateU );
   WALBERLA_LOG_INFO_ON_ROOT( "  - l2 error p:    " << std::scientific << avgl2ErrorConvergenceRateP );
   WALBERLA_LOG_INFO_ON_ROOT( "  - l2 residual p: " << std::scientific << avgl2ResidualConvergenceRateP );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

void setup( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;
   LIKWID_MARKER_REGISTER( "FMG" );
   LIKWID_MARKER_REGISTER( "VCYCLE" );

#ifdef HHG_BUILD_WITH_PETSC
   PETScManager petscManager;
#endif

   WALBERLA_LOG_INFO_ON_ROOT( "///////////////////////" );
   WALBERLA_LOG_INFO_ON_ROOT( "// Multigrid Studies //" );
   WALBERLA_LOG_INFO_ON_ROOT( "///////////////////////" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./MultigridStudies.prm";
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   ////////////////
   // Parameters //
   ////////////////

   const std::string equation                        = mainConf.getParameter< std::string >( "equation" );
   const uint_t      dim                             = mainConf.getParameter< uint_t >( "dim" );
   const uint_t      numProcesses                    = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
   const uint_t      numFacesPerSide                 = mainConf.getParameter< uint_t >( "numFacesPerSide" );
   const std::string discretization                  = mainConf.getParameter< std::string >( "discretization" );
   const uint_t      numCycles                       = mainConf.getParameter< uint_t >( "numCycles" );
   const std::string cycleTypeString                 = mainConf.getParameter< std::string >( "cycleType" );
   const uint_t      fmgInnerCycles                  = mainConf.getParameter< uint_t >( "fmgInnerCycles" );
   const real_t      L2residualTolerance             = mainConf.getParameter< real_t >( "L2residualTolerance" );
   const real_t      sorRelax                        = mainConf.getParameter< real_t >( "sorRelax" );
   const bool        symmGSVelocity                  = mainConf.getParameter< bool >( "symmGSVelocity" );
   const bool        symmGSPressure                  = mainConf.getParameter< bool >( "symmGSPressure" );
   const uint_t      numGSVelocity                   = mainConf.getParameter< uint_t >( "numGSVelocity" );
   const uint_t      numGSPressure                   = mainConf.getParameter< uint_t >( "numGSPressure" );
   const uint_t      preSmoothingSteps               = mainConf.getParameter< uint_t >( "preSmoothingSteps" );
   const uint_t      postSmoothingSteps              = mainConf.getParameter< uint_t >( "postSmoothingSteps" );
   const uint_t      smoothingIncrement              = mainConf.getParameter< uint_t >( "smoothingIncrement" );
   const bool        projectPressureAfterRestriction = mainConf.getParameter< bool >( "projectPressureAfterRestriction" );
   const uint_t      minLevel                        = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel                        = ( discretization == "P1" ? mainConf.getParameter< uint_t >( "maxLevel" ) :
                                                      mainConf.getParameter< uint_t >( "maxLevel" ) - 1 );
   const bool        calculateDiscretizationError    = mainConf.getParameter< bool >( "calculateDiscretizationError" );
   const uint_t      coarseGridMaxIterations         = mainConf.getParameter< uint_t >( "coarseGridMaxIterations" );
   const real_t      coarseGridResidualTolerance     = mainConf.getParameter< real_t >( "coarseGridResidualTolerance" );
   const bool        outputVTK                       = mainConf.getParameter< bool >( "outputVTK" );
   const bool        outputTiming                    = mainConf.getParameter< bool >( "outputTiming" );
   const bool        outputTimingJSON                = mainConf.getParameter< bool >( "outputTimingJSON" );
   const std::string outputTimingJSONFile            = mainConf.getParameter< std::string >( "outputTimingJSONFile" );
   const bool        outputSQL                       = mainConf.getParameter< bool >( "outputSQL" );
   const std::string outputSQLFile                   = mainConf.getParameter< std::string >( "outputSQLFile" );
   const std::string sqlTag                          = mainConf.getParameter< std::string >( "sqlTag", "default" );
   const uint_t      skipCyclesForAvgConvRate        = mainConf.getParameter< uint_t >( "skipCyclesForAvgConvRate" );
   const std::string meshLayout                      = mainConf.getParameter< std::string >( "meshLayout" );
   const bool        symmetricCuboidMesh             = mainConf.getParameter< bool >( "symmetricCuboidMesh" );
   const uint_t      cyclesBeforeDC                  = mainConf.getParameter< uint_t >( "cyclesBeforeDC" );
   const uint_t      postDCPreSmoothingSteps         = mainConf.getParameter< uint_t >( "postDCPreSmoothingSteps" );
   const uint_t      postDCPostSmoothingSteps        = mainConf.getParameter< uint_t >( "postDCPostSmoothingSteps" );
   const uint_t      postDCSmoothingIncrement        = mainConf.getParameter< uint_t >( "postDCSmoothingIncrement" );

   const bool        meshSphericalShell              = mainConf.getParameter< bool >( "meshSphericalShell" );
   const uint_t      shellNTan                       = mainConf.getParameter< uint_t >( "shellNTan" );
   const uint_t      shellNRad                       = mainConf.getParameter< uint_t >( "shellNRad" );
   const real_t      shellRMin                       = mainConf.getParameter< real_t >( "shellRMin" );
   const real_t      shellRMax                       = mainConf.getParameter< real_t >( "shellRMax" );

   // parameter checks
   WALBERLA_CHECK( equation == "stokes" || equation == "poisson" );
   WALBERLA_CHECK( discretization == "P1" || discretization == "P2" );
   WALBERLA_CHECK( cycleTypeString == "V" || cycleTypeString == "W" );

   const CycleType cycleType = ( cycleTypeString == "V" ? CycleType::VCYCLE : CycleType::WCYCLE );

   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - equation:                                " << equation );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num processes:                           " << numProcesses );
   WALBERLA_LOG_INFO_ON_ROOT( "  - discretization:                          " << discretization );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num cycles:                              " << numCycles );
   WALBERLA_LOG_INFO_ON_ROOT( "  - cycle type:                              " << cycleTypeString );
   WALBERLA_LOG_INFO_ON_ROOT(
       "  - full multigrid:                          "
       << ( fmgInnerCycles == 0 ? "no" : "yes, inner cycles per level: " + std::to_string( fmgInnerCycles ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual tolerance:                   " << L2residualTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "  - SOR relax:                               " << sorRelax );
   WALBERLA_LOG_INFO_ON_ROOT( "  - Uzawa velocity smoother:                 " << ( symmGSVelocity ? "symmetric" : "forward" )
                                                                              << " GS, " << numGSVelocity << " iterations" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - Uzawa pressure smoother:                 " << ( symmGSPressure ? "symmetric" : "forward" )
                                                                              << " SOR, " << numGSPressure << " iterations" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - pre- / post- / incr-smoothing:           " << preSmoothingSteps << " / " << postSmoothingSteps
                                                                              << " / " << smoothingIncrement );
   WALBERLA_LOG_INFO_ON_ROOT( "  - min / max level:                         " << minLevel << " / " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT(
       "  - project pressure after restriction:      " << ( projectPressureAfterRestriction ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT(
       "  - calculate discretization error:          " << ( calculateDiscretizationError ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - coarse grid max itertions (stokes only): " << coarseGridMaxIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "  - coarse grid residual tol  (stokes only): " << coarseGridResidualTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output VTK:                              " << ( outputVTK ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing:                           " << ( outputTiming ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing JSON:                      " << ( outputTimingJSON ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing JSON file:                 " << outputTimingJSONFile );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output SQL:                              " << ( outputSQL ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output SQL file:                         " << outputSQLFile );
   WALBERLA_LOG_INFO_ON_ROOT( "  - SQL tag:                                 " << sqlTag );
   WALBERLA_LOG_INFO_ON_ROOT( "  - skip cycles for avg conv rate:           " << skipCyclesForAvgConvRate );
   if ( meshSphericalShell )
   {
     WALBERLA_LOG_INFO_ON_ROOT( "  - sphericalShell:                          " << "yes" );
     WALBERLA_LOG_INFO_ON_ROOT( "  - nTan:                                    " << shellNTan );
     WALBERLA_LOG_INFO_ON_ROOT( "  - nRad:                                    " << shellNRad );
     WALBERLA_LOG_INFO_ON_ROOT( "  - rMin:                                    " << shellRMin );
     WALBERLA_LOG_INFO_ON_ROOT( "  - rMax:                                    " << shellRMax );
   }
   else
   {
     WALBERLA_LOG_INFO_ON_ROOT( "  - dim:                                     " << dim );
     WALBERLA_LOG_INFO_ON_ROOT( "  - num faces per side:                      " << numFacesPerSide );
     WALBERLA_LOG_INFO_ON_ROOT( "  - mesh layout:                             " << meshLayout );
     WALBERLA_LOG_INFO_ON_ROOT( "  - symmetric cuboid mesh:                   " << symmetricCuboidMesh );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "  - cycles before DC:                        "
                              << ( discretization == "P1" ? std::to_string( cyclesBeforeDC ) : "disabled" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - DC pre- / post- / incr-smoothing:        "
                              << ( discretization == "P1" ? std::to_string( postDCPreSmoothingSteps ) + " / " +
                                                                std::to_string( postDCPostSmoothingSteps ) + " / " +
                                                                std::to_string( postDCSmoothingIncrement ) :
                                                            "disabled" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   /////////
   // SQL //
   /////////

   std::map< std::string, walberla::int64_t >          sqlIntegerProperties;
   std::map< std::string, double >                     sqlRealProperties;
   std::map< std::string, std::string >                sqlStringProperties;
   std::map< uint_t, std::map< std::string, double > > sqlRealPropertiesMG;

   sqlStringProperties["tag"]      = sqlTag;
   sqlStringProperties["equation"] = equation;

   sqlIntegerProperties["num_processes"]               = int64_c( numProcesses );
   sqlIntegerProperties["dim"]                         = int64_c( dim );
   sqlIntegerProperties["num_faces_per_side"]          = int64_c( numFacesPerSide );
   sqlStringProperties["discretization"]               = discretization;
   sqlIntegerProperties["num_cycles"]                  = int64_c( numCycles );
   sqlStringProperties["cycle_type"]                   = cycleTypeString;
   sqlIntegerProperties["fmgInnerCycles"]              = int64_c( fmgInnerCycles );
   sqlRealProperties["sor_relax"]                      = sorRelax;
   sqlIntegerProperties["pre_smoothing"]               = int64_c( preSmoothingSteps );
   sqlIntegerProperties["post_smoothing"]              = int64_c( postSmoothingSteps );
   sqlIntegerProperties["incr_smoothing"]              = int64_c( smoothingIncrement );
   sqlIntegerProperties["min_level"]                   = int64_c( minLevel );
   sqlIntegerProperties["max_level"]                   = int64_c( maxLevel );
   sqlIntegerProperties["coarse_grid_max_iterations"]  = int64_c( coarseGridMaxIterations );
   sqlRealProperties["coarse_grid_residual_tolerance"] = coarseGridResidualTolerance;
   sqlIntegerProperties["project_after_restriction"]   = int64_c( projectPressureAfterRestriction );
   sqlIntegerProperties["symm_gs_velocity"]            = int64_c( symmGSVelocity );
   sqlIntegerProperties["symm_gs_pressure"]            = int64_c( symmGSPressure );
   sqlIntegerProperties["num_gs_velocity"]             = int64_c( numGSVelocity );
   sqlIntegerProperties["num_gs_pressure"]             = int64_c( numGSPressure );

   sqlIntegerProperties["cycles_before_dc"]  = int64_c( cyclesBeforeDC );
   sqlIntegerProperties["dc_pre_smoothing"]  = int64_c( preSmoothingSteps );
   sqlIntegerProperties["dc_post_smoothing"] = int64_c( postSmoothingSteps );
   sqlIntegerProperties["dc_incr_smoothing"] = int64_c( smoothingIncrement );

   ////////////
   // Domain //
   ////////////

   Point2D leftBottom( {0, 0} );
   Point3D leftBottom3D( {0, 0, 0} );
   if ( equation == "stokes" && ( NEUMANN_PROBLEM || COLLIDING_FLOW ) )
   {
      leftBottom   = Point2D( {-1, -1} );
      leftBottom3D = Point3D( {-1, -1, -1} );
   }

   std::shared_ptr< PrimitiveStorage > storage;
   {
     if ( meshSphericalShell )
     {
       auto meshInfo = MeshInfo::meshSphericalShell( shellNTan, shellNRad, shellRMin, shellRMax );

       SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
       setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

       storage = std::make_shared< PrimitiveStorage >( setupStorage );

       sqlIntegerProperties["num_macro_vertices"] = int64_c( setupStorage.getNumberOfVertices() );
       sqlIntegerProperties["num_macro_edges"]    = int64_c( setupStorage.getNumberOfEdges() );
       sqlIntegerProperties["num_macro_faces"]    = int64_c( setupStorage.getNumberOfFaces() );
       sqlIntegerProperties["num_macro_cells"]    = int64_c( setupStorage.getNumberOfCells() );

       exactU = shellExactU;
       exactV = shellExactV;
       exactW = shellExactW;

       exactP = shellExactP;

       bcU = shellExactU;
       bcV = shellExactV;
       bcW = shellExactW;

       rhsU = shellRhsU;
       rhsV = shellRhsV;
       rhsW = shellRhsW;
     }
     else
     {
       MeshInfo::meshFlavour meshFlavour;
       if ( meshLayout == "CRISS" )
         meshFlavour = MeshInfo::CRISS;
       else if ( meshLayout == "CRISSCROSS" )
         meshFlavour = MeshInfo::CRISSCROSS;
       else
       WALBERLA_ABORT( "Invalid mesh layout." );

       auto meshInfo = MeshInfo::meshRectangle( leftBottom, Point2D( { 1, 1 } ), meshFlavour, numFacesPerSide, numFacesPerSide );
       if ( dim == 3 )
       {
         if ( symmetricCuboidMesh )
         {
           meshInfo = MeshInfo::meshSymmetricCuboid(
           leftBottom3D, Point3D( { 1, 1, 1 } ), numFacesPerSide, numFacesPerSide, numFacesPerSide );
         } else
         {
           meshInfo =
           MeshInfo::meshCuboid( leftBottom3D, Point3D( { 1, 1, 1 } ), numFacesPerSide, numFacesPerSide, numFacesPerSide );
         }
       }
       SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
       setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

#if NEUMANN_PROBLEM
       auto topBoundary = []( const Point3D& x ) {
          const real_t eps = 1e-8;
          return std::abs( x[1] - 1 ) < eps;
       };
       auto bottomBoundary = []( const Point3D& x ) {
          const real_t eps = 1e-8;
          return std::abs( x[1] + 1 ) < eps;
       };
       auto leftBoundary = []( const Point3D& x ) {
          const real_t eps = 1e-8;
          return std::abs( x[0] + 1 ) < eps;
       };
       setupStorage.setMeshBoundaryFlagsOnBoundary( 2, 0, true );
       setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, topBoundary, true );
       setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, bottomBoundary, true );
       setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, leftBoundary, true );
#endif

#if CONSTANTA_POISSON
       if ( dim == 2 )
       {
         exact = exactConstanta2D;
         rhs = rhsConstanta2D;
       } else
       {
         exact = exactConstanta3D;
         rhs = rhsConstanta3D;
       }
#endif
       storage = std::make_shared< PrimitiveStorage >( setupStorage );

       sqlIntegerProperties["num_macro_vertices"] = int64_c( setupStorage.getNumberOfVertices() );
       sqlIntegerProperties["num_macro_edges"]    = int64_c( setupStorage.getNumberOfEdges() );
       sqlIntegerProperties["num_macro_faces"]    = int64_c( setupStorage.getNumberOfFaces() );
       sqlIntegerProperties["num_macro_cells"]    = int64_c( setupStorage.getNumberOfCells() );
     }

   }

   if ( outputVTK )
   {
      writeDomainPartitioningVTK( storage, "vtk", "Domain" );
   }

   auto globalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( globalInfo );

   if ( equation == "poisson" )
   {
      if ( discretization == "P1" )
      {
         MultigridLaplace< P1Function< real_t >,
                           P1ConstantLaplaceOperator,
                           P1ConstantMassOperator,
                           P1toP1LinearRestriction,
                           P1toP1LinearProlongation,
                           P1toP1QuadraticProlongation >( storage,
                                                          minLevel,
                                                          maxLevel,
                                                          numCycles,
                                                          cycleType,
                                                          fmgInnerCycles,
                                                          L2residualTolerance,
                                                          sorRelax,
                                                          preSmoothingSteps,
                                                          postSmoothingSteps,
                                                          outputVTK,
                                                          skipCyclesForAvgConvRate,
                                                          calculateDiscretizationError,
                                                          sqlIntegerProperties,
                                                          sqlRealProperties,
                                                          sqlStringProperties,
                                                          sqlRealPropertiesMG );
      }
      else if ( discretization == "P2" )
      {
         MultigridLaplace< P2Function< real_t >,
                           P2ConstantLaplaceOperator,
                           P2ConstantMassOperator,
                           P2toP2QuadraticRestriction,
                           P2toP2QuadraticProlongation,
                           P2toP2QuadraticProlongation >( storage,
                                                          minLevel,
                                                          maxLevel,
                                                          numCycles,
                                                          cycleType,
                                                          fmgInnerCycles,
                                                          L2residualTolerance,
                                                          sorRelax,
                                                          preSmoothingSteps,
                                                          postSmoothingSteps,
                                                          outputVTK,
                                                          skipCyclesForAvgConvRate,
                                                          calculateDiscretizationError,
                                                          sqlIntegerProperties,
                                                          sqlRealProperties,
                                                          sqlStringProperties,
                                                          sqlRealPropertiesMG );
      }
   }
   else if ( equation == "stokes" )
   {
      if ( discretization == "P1" )
      {
         MultigridStokes< P1StokesFunction< real_t >,
                          P1StokesFunction< int >,
                          P1StokesOperator,
                          P1ConstantMassOperator,
                          P1P1StokesToP1P1StokesRestriction,
                          P1P1StokesToP1P1StokesProlongation,
                          P1P1StokesToP1P1StokesProlongation >( storage,
                                                                minLevel,
                                                                maxLevel,
                                                                numCycles,
                                                                cycleType,
                                                                fmgInnerCycles,
                                                                L2residualTolerance,
                                                                sorRelax,
                                                                symmGSVelocity,
                                                                numGSVelocity,
                                                                symmGSPressure,
                                                                numGSPressure,
                                                                preSmoothingSteps,
                                                                postSmoothingSteps,
                                                                smoothingIncrement,
                                                                projectPressureAfterRestriction,
                                                                coarseGridMaxIterations,
                                                                coarseGridResidualTolerance,
                                                                outputVTK,
                                                                skipCyclesForAvgConvRate,
                                                                calculateDiscretizationError,
                                                                cyclesBeforeDC,
                                                                postDCPreSmoothingSteps,
                                                                postDCPostSmoothingSteps,
                                                                postDCSmoothingIncrement,
                                                                sqlIntegerProperties,
                                                                sqlRealProperties,
                                                                sqlStringProperties,
                                                                sqlRealPropertiesMG );
      }
      else if ( discretization == "P2" )
      {
         MultigridStokes< P2P1TaylorHoodFunction< real_t >,
                          P2P1TaylorHoodFunction< int >,
                          P2P1TaylorHoodStokesOperator,
                          P2ConstantMassOperator,
                          P2P1StokesToP2P1StokesRestriction,
                          P2P1StokesToP2P1StokesProlongation,
                          P2P1StokesToP2P1StokesProlongation >( storage,
                                                                minLevel,
                                                                maxLevel,
                                                                numCycles,
                                                                cycleType,
                                                                fmgInnerCycles,
                                                                L2residualTolerance,
                                                                sorRelax,
                                                                symmGSVelocity,
                                                                numGSVelocity,
                                                                symmGSPressure,
                                                                numGSPressure,
                                                                preSmoothingSteps,
                                                                postSmoothingSteps,
                                                                smoothingIncrement,
                                                                projectPressureAfterRestriction,
                                                                coarseGridMaxIterations,
                                                                coarseGridResidualTolerance,
                                                                outputVTK,
                                                                skipCyclesForAvgConvRate,
                                                                calculateDiscretizationError,
                                                                0,
                                                                preSmoothingSteps,
                                                                postSmoothingSteps,
                                                                smoothingIncrement,
                                                                sqlIntegerProperties,
                                                                sqlRealProperties,
                                                                sqlStringProperties,
                                                                sqlRealPropertiesMG );
      }
   }

   auto tt = storage->getTimingTree()->getReduced().getCopyWithRemainder();
   if ( outputTiming )
   {
      WALBERLA_LOG_INFO_ON_ROOT( tt );
   }

   WALBERLA_ROOT_SECTION()
   {
      if ( outputTimingJSON )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Writing JSON timing to " << outputTimingJSONFile )
         nlohmann::json ttJson;
         walberla::timing::to_json( ttJson, tt );
         std::ofstream jsonOutput;
         jsonOutput.open( outputTimingJSONFile );
         jsonOutput << ttJson.dump( 4 );
         jsonOutput.close();
         WALBERLA_LOG_INFO_ON_ROOT( "Done writing JSON timing." )
      }
   }

   if ( outputSQL )
   {
      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Writing SQL database to " << outputSQLFile )
         const std::string                  dbFile = outputSQLFile;
         walberla::sqlite::SQLiteDB db( dbFile );
         sqlIntegerProperties["conv_table_for_run"] = -1;
         auto runId                                 = db.storeRun( sqlIntegerProperties, sqlStringProperties, sqlRealProperties );
         for ( uint_t cycle = 0; cycle <= numCycles; cycle++ )
         {
            if ( sqlRealPropertiesMG.count( cycle ) > 0 )
            {
               std::map< std::string, int64_t > runIdMap;
               runIdMap["conv_table_for_run"] = int64_c( runId );
               runIdMap["cycle"]              = int64_c( cycle );
               db.storeRun( runIdMap, std::map< std::string, std::string >(), sqlRealPropertiesMG[cycle] );
            }
         }
         WALBERLA_LOG_INFO_ON_ROOT( "Done writing SQL database." )
      }
   }

   LIKWID_MARKER_CLOSE;
}

} // namespace hyteg

int main( int argc, char** argv )
{
   hyteg::setup( argc, argv );
}
