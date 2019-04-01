
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/timing/TimingJSON.h"
#include "core/math/Constants.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"
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

#include "postprocessing/sqlite/SQLite.h"

namespace hhg {

using walberla::int64_c;
using walberla::math::M_PI;

#define NEUMANN_PROBLEM 0
#define COLLIDING_FLOW 0

#if 0
std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x )
{
   return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
};

std::function< real_t( const hhg::Point3D& ) > rhs = []( const hhg::Point3D& x )
{
   return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
};
#else
std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };

std::function< real_t( const hhg::Point3D& ) > rhs    = []( const hhg::Point3D& ) { return 0; };
#endif

#if NEUMANN_PROBLEM
std::function< real_t( const hhg::Point3D& ) > bcU = []( const hhg::Point3D& x ) {
   if ( std::abs( x[0] + 1 ) < 1e-8 )
   {
      return 1 - x[1] * x[1];
   }
   else
   {
      return 0.0;
   }
};

std::function< real_t( const hhg::Point3D& ) > exactU = []( const hhg::Point3D& x ) { return 1 - x[1] * x[1]; };

std::function< real_t( const hhg::Point3D& ) > exactV = []( const hhg::Point3D& ) { return 0.0; };
std::function< real_t( const hhg::Point3D& ) > exactP = []( const hhg::Point3D& x ) { return -2 * x[0]; };
std::function< real_t( const hhg::Point3D& ) > rhsU   = []( const hhg::Point3D& ) { return 0; };
std::function< real_t( const hhg::Point3D& ) > rhsV   = []( const hhg::Point3D& x ) { return 0; };

#else

#if COLLIDING_FLOW
std::function< real_t( const hhg::Point3D& ) > exactU = []( const hhg::Point3D& x ) { return 20 * x[0] * std::pow( x[1], 3.0 ); };
std::function< real_t( const hhg::Point3D& ) > bcU    = []( const hhg::Point3D& x ) { return 20 * x[0] * std::pow( x[1], 3.0 ); };
std::function< real_t( const hhg::Point3D& ) > exactV = []( const hhg::Point3D& x ) {
   return 5 * std::pow( x[0], 4.0 ) - 5 * std::pow( x[1], 4.0 );
};
std::function< real_t( const hhg::Point3D& ) > exactP = []( const hhg::Point3D& x ) {
   return 60 * std::pow( x[0], 2.0 ) * x[1] - 20 * std::pow( x[1], 3.0 );
};
std::function< real_t( const hhg::Point3D& ) > rhsU = []( const hhg::Point3D& ) { return 0; };
std::function< real_t( const hhg::Point3D& ) > rhsV = []( const hhg::Point3D& x ) { return 0; };

#else
std::function< real_t( const hhg::Point3D& ) > exactU = []( const hhg::Point3D& x ) {
   return std::sin( 2 * M_PI * x[0] ) * std::cos( M_PI * x[1] );
};
std::function< real_t( const hhg::Point3D& ) > bcU = []( const hhg::Point3D& x ) {
   return std::sin( 2 * M_PI * x[0] ) * std::cos( M_PI * x[1] );
};
std::function< real_t( const hhg::Point3D& ) > exactV = []( const hhg::Point3D& x ) {
   return -2.0 * std::cos( 2 * M_PI * x[0] ) * std::sin( M_PI * x[1] );
};
std::function< real_t( const hhg::Point3D& ) > exactP = []( const hhg::Point3D& x ) {
   return 2.5 * M_PI * std::cos( 2 * M_PI * x[0] ) * std::cos( M_PI * x[1] );
};
std::function< real_t( const hhg::Point3D& ) > rhsU = []( const hhg::Point3D& ) { return 0; };
std::function< real_t( const hhg::Point3D& ) > rhsV = []( const hhg::Point3D& x ) {
   return -12.5 * M_PI * M_PI * std::cos( 2 * M_PI * x[0] ) * std::sin( M_PI * x[1] );
};
#endif
#endif

template < typename Function, typename LaplaceOperator, typename MassOperator >
void calculateErrorAndResidual( const uint_t&          level,
                                const LaplaceOperator& A,
                                const MassOperator&    M,
                                const Function&        u,
                                const Function&        f,
                                const Function&        uExact,
                                const Function&        error,
                                const Function&        residual,
                                const Function&        tmp,
                                real_t&                l2Error,
                                real_t&                L2Error,
                                real_t&                l2Residual,
                                real_t&                L2Residual )
{
   error.assign( {1.0, -1.0}, {uExact, u}, level, All );

   tmp.interpolate( real_c( 0 ), level, All );
   A.apply( u, tmp, level, Inner );
   residual.assign( {1.0, -1.0}, {f, tmp}, level, All );

   M.apply( error, tmp, level, Inner );
   l2Error = std::sqrt( error.dotGlobal( error, level, Inner ) );
   L2Error = std::sqrt( error.dotGlobal( tmp, level, Inner ) );
   M.apply( residual, tmp, level, Inner );
   l2Residual = std::sqrt( residual.dotGlobal( residual, level, Inner ) );
   L2Residual = std::sqrt( residual.dotGlobal( tmp, level, Inner ) );
   L2Error    = l2Error;
   L2Residual = l2Residual;
}

template < typename Function, typename StokesOperator >
void calculateErrorAndResidualStokes( const uint_t&         level,
                                      const StokesOperator& A,
                                      const Function&       u,
                                      const Function&       f,
                                      const Function&       uExact,
                                      const Function&       error,
                                      const Function&       residual,
                                      const Function&       tmp,
                                      real_t&               l2ErrorU,
                                      real_t&               l2ErrorV,
                                      real_t&               l2ErrorP,
                                      real_t&               l2ResidualU,
                                      real_t&               l2ResidualV,
                                      real_t&               l2ResidualP )
{
   error.assign( {1.0, -1.0}, {uExact, u}, level, All );

   tmp.interpolate( real_c( 0 ), level, All );
   A.apply( u, tmp, level, Inner | NeumannBoundary );
   residual.assign( {1.0, -1.0}, {f, tmp}, level, All );

   vertexdof::projectMean( error.p, level );

   l2ErrorU    = std::sqrt( error.u.dotGlobal( error.u, level, Inner | NeumannBoundary ) );
   l2ErrorV    = std::sqrt( error.v.dotGlobal( error.v, level, Inner | NeumannBoundary ) );
   l2ErrorP    = std::sqrt( error.p.dotGlobal( error.p, level, Inner | NeumannBoundary ) );
   l2ResidualU = std::sqrt( residual.u.dotGlobal( residual.u, level, Inner | NeumannBoundary ) );
   l2ResidualV = std::sqrt( residual.v.dotGlobal( residual.v, level, Inner | NeumannBoundary ) );
   l2ResidualP = std::sqrt( residual.p.dotGlobal( residual.p, level, Inner | NeumannBoundary ) );
}

template < typename Function, typename LaplaceOperator, typename MassOperator >
void calculateDiscretizationError( const std::shared_ptr< PrimitiveStorage >& storage,
                                   const uint_t&                              level,
                                   real_t&                                    l2DiscretizationError )
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

   auto solver = std::make_shared< CGSolver< LaplaceOperator > >( storage, level, level );
   solver->solve( A, u, f, level );

   real_t L2Error;
   real_t l2Residual;
   real_t L2Residual;
   calculateErrorAndResidual(
       level, A, M, u, f, uExact, error, residual, tmp, l2DiscretizationError, L2Error, l2Residual, L2Residual );
}

template < typename StokesFunction, typename StokesOperator, typename MassOperator >
void calculateDiscretizationErrorStokes( const std::shared_ptr< PrimitiveStorage >& storage,
                                         const uint_t&                              level,
                                         real_t&                                    l2DiscretizationErrorU,
                                         real_t&                                    l2DiscretizationErrorV,
                                         real_t&                                    l2DiscretizationErrorP )
{
   StokesFunction u( "u", storage, level, level );
   StokesFunction f( "f", storage, level, level );

   StokesFunction uExact( "uExact", storage, level, level );
   StokesFunction residual( "residual", storage, level, level );
   StokesFunction error( "error", storage, level, level );
   StokesFunction tmp( "tmp", storage, level, level );

   StokesOperator A( storage, level, level );
   MassOperator   M( storage, level, level );

   u.u.interpolate( bcU, level, DirichletBoundary );
   u.v.interpolate( exactV, level, DirichletBoundary );

   uExact.u.interpolate( exactU, level, All );
   uExact.v.interpolate( exactV, level, All );
   uExact.p.interpolate( exactP, level, All );

   tmp.u.interpolate( rhsU, level, All );
   tmp.v.interpolate( rhsV, level, All );
   M.apply( tmp.u, f.u, level, All );
   M.apply( tmp.v, f.v, level, All );

#ifdef HHG_BUILD_WITH_PETSC
   auto solver = std::make_shared< PETScLUSolver< StokesOperator > >( storage, level );
#else
   auto cgVelocity =
       std::make_shared< CGSolver< typename StokesOperator::VelocityOperator_T > >( storage, level, level, 0, 1e-14 );
   auto preconditioner = std::make_shared< StokesBlockDiagonalPreconditioner< StokesOperator, P1LumpedInvMassOperator > >(
       storage, level, level, 1, cgVelocity );
   auto solver = std::make_shared< MinResSolver< StokesOperator > >(
       storage, level, level, std::numeric_limits< uint_t >::max(), 1e-16, preconditioner );
#endif
   solver->solve( A, u, f, level );

   vertexdof::projectMean( u.p, level );

   real_t l2ResidualU;
   real_t l2ResidualV;
   real_t l2ResidualP;

   calculateErrorAndResidualStokes( level,
                                    A,
                                    u,
                                    f,
                                    uExact,
                                    error,
                                    residual,
                                    tmp,
                                    l2DiscretizationErrorU,
                                    l2DiscretizationErrorV,
                                    l2DiscretizationErrorP,
                                    l2ResidualU,
                                    l2ResidualV,
                                    l2ResidualP );
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

   real_t l2Error;
   real_t L2Error;
   real_t l2Residual;
   real_t L2Residual;

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

   if ( calcDiscretizationError )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "l2 discretization error per level:" );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         real_t discretizationError;
         calculateDiscretizationError< Function, LaplaceOperator, MassOperator >( storage, level, discretizationError );
         WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::scientific << discretizationError );
         sqlRealProperties["l2_discr_error_level_" + std::to_string( level )] = discretizationError;
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
   calculateErrorAndResidual( maxLevel, A, M, u, f, uExact, error, residual, tmp, l2Error, L2Error, l2Residual, L2Residual );
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
       " After cycle... ||     l2 error |     L2 error | L2 error reduction ||  l2 residual |  L2 residual | L2 residual reduction || time cycle [s] | time error calculation [s] | time VTK [s] |" );
   WALBERLA_LOG_INFO_ON_ROOT(
       " ---------------++--------------+--------------+--------------------++--------------+--------------+-----------------------++----------------+----------------------------+--------------|" );
   WALBERLA_LOG_INFO_ON_ROOT( "        initial || " << std::scientific << l2Error << " | " << L2Error << " | "
                                                    << "               --- || " << l2Residual << " | " << L2Residual
                                                    << " |                   --- ||            --- | " << std::fixed
                                                    << std::setprecision( 2 ) << std::setw( 26 ) << timeError << " | "
                                                    << std::setw( 12 ) << timeVTK << " |" );

   real_t avgL2ErrorConvergenceRate    = 0;
   real_t avgL2ResidualConvergenceRate = 0;
   real_t avgl2ErrorConvergenceRate    = 0;
   real_t avgl2ResidualConvergenceRate = 0;

   real_t L2ErrorReduction    = 0;
   real_t L2ResidualReduction = 0;
   real_t l2ErrorReduction    = 0;
   real_t l2ResidualReduction = 0;

   sqlRealPropertiesMG[0]["capital_L2_error"]              = L2Error;
   sqlRealPropertiesMG[0]["capital_L2_error_reduction"]    = L2ErrorReduction;
   sqlRealPropertiesMG[0]["capital_L2_residual"]           = L2Residual;
   sqlRealPropertiesMG[0]["capital_L2_residual_reduction"] = L2ResidualReduction;

   sqlRealPropertiesMG[0]["lowercase_l2_error"]              = l2Error;
   sqlRealPropertiesMG[0]["lowercase_l2_error_reduction"]    = l2ErrorReduction;
   sqlRealPropertiesMG[0]["lowercase_l2_residual"]           = l2Residual;
   sqlRealPropertiesMG[0]["lowercase_l2_residual_reduction"] = l2ResidualReduction;

   ///////////
   // Solve //
   ///////////

   auto smoother         = std::make_shared< SORSmoother< LaplaceOperator > >( sorRelax );
   auto coarseGridSolver = std::make_shared< CGSolver< LaplaceOperator > >( storage, minLevel, minLevel );

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
      real_t _l2Error, _L2Error, _l2Residual, _L2Residual;
      calculateErrorAndResidual(
          currentLevel, A, M, u, f, uExact, error, residual, tmp, _l2Error, _L2Error, _l2Residual, _L2Residual );
      sqlRealProperties["fmg_l2_error_level_" + std::to_string( currentLevel )] = _l2Error;
      WALBERLA_LOG_INFO_ON_ROOT( "    fmg level " << currentLevel << ": l2 error: " << std::scientific << _l2Error );
   };

   FullMultigridSolver< LaplaceOperator > fullMultigridSolver(
       storage, multigridSolver, fmgProlongation, minLevel, maxLevel, fmgInnerCycles, postCycle );

   for ( uint_t cycle = 1; cycle <= numCycles; cycle++ )
   {
      const real_t lastL2Error    = L2Error;
      const real_t lastL2Residual = L2Residual;

      const real_t lastl2Error    = l2Error;
      const real_t lastl2Residual = l2Residual;

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

      timer.reset();
      calculateErrorAndResidual( maxLevel, A, M, u, f, uExact, error, residual, tmp, l2Error, L2Error, l2Residual, L2Residual );
      timer.end();
      timeError = timer.last();

      timer.reset();
      if ( outputVTK )
      {
         vtkOutput.write( maxLevel, cycle );
      }
      timer.end();
      timeVTK = timer.last();

      L2ErrorReduction    = L2Error / lastL2Error;
      L2ResidualReduction = L2Residual / lastL2Residual;
      l2ErrorReduction    = l2Error / lastl2Error;
      l2ResidualReduction = l2Residual / lastl2Residual;

      WALBERLA_LOG_INFO_ON_ROOT( std::setw( 15 ) << cycle << " || " << std::scientific << l2Error << " | " << L2Error << " | "
                                                 << "      " << L2ErrorReduction << " || " << l2Residual << " | " << L2Residual
                                                 << " |          " << L2ResidualReduction << " || " << std::fixed
                                                 << std::setprecision( 2 ) << std::setw( 14 ) << timeCycle << " | "
                                                 << std::setw( 26 ) << timeError << " | " << std::setw( 12 ) << timeVTK << " |" );

      if ( cycle > skipCyclesForAvgConvRate )
      {
         avgL2ErrorConvergenceRate += L2ErrorReduction;
         avgL2ResidualConvergenceRate += L2ResidualReduction;
         avgl2ErrorConvergenceRate += l2ErrorReduction;
         avgl2ResidualConvergenceRate += l2ResidualReduction;
      }

      sqlRealPropertiesMG[cycle]["capital_L2_error"]              = L2Error;
      sqlRealPropertiesMG[cycle]["capital_L2_error_reduction"]    = L2ErrorReduction;
      sqlRealPropertiesMG[cycle]["capital_L2_residual"]           = L2Residual;
      sqlRealPropertiesMG[cycle]["capital_L2_residual_reduction"] = L2ResidualReduction;

      sqlRealPropertiesMG[cycle]["lowercase_l2_error"]              = l2Error;
      sqlRealPropertiesMG[cycle]["lowercase_l2_error_reduction"]    = l2ErrorReduction;
      sqlRealPropertiesMG[cycle]["lowercase_l2_residual"]           = l2Residual;
      sqlRealPropertiesMG[cycle]["lowercase_l2_residual_reduction"] = l2ResidualReduction;

      if ( L2Residual < L2residualTolerance )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "L2 residual dropped below tolerance." )
         break;
      }
   }

   avgL2ErrorConvergenceRate /= real_c( numCycles - skipCyclesForAvgConvRate );
   avgL2ResidualConvergenceRate /= real_c( numCycles - skipCyclesForAvgConvRate );

   avgl2ErrorConvergenceRate /= real_c( numCycles - skipCyclesForAvgConvRate );
   avgl2ResidualConvergenceRate /= real_c( numCycles - skipCyclesForAvgConvRate );

   sqlRealProperties["avg_capital_L2_error_conv_rate"]    = avgL2ErrorConvergenceRate;
   sqlRealProperties["avg_capital_L2_residual_conv_rate"] = avgL2ResidualConvergenceRate;

   sqlRealProperties["avg_lowercase_l2_error_conv_rate"]    = avgl2ErrorConvergenceRate;
   sqlRealProperties["avg_lowercase_l2_residual_conv_rate"] = avgl2ResidualConvergenceRate;

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Average convergence rates:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 error:    " << std::scientific << avgL2ErrorConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual: " << std::scientific << avgL2ResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - l2 error:    " << std::scientific << avgl2ErrorConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - l2 residual: " << std::scientific << avgl2ResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
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
                      const uint_t&                                        preSmoothingSteps,
                      const uint_t&                                        postSmoothingSteps,
                      const uint_t&                                        smoothingIncrement,
                      const bool&                                          projectPressureAfterRestriction,
                      const uint_t&                                        coarseGridMaxIterations,
                      const bool&                                          outputVTK,
                      const uint_t&                                        skipCyclesForAvgConvRate,
                      const bool&                                          calcDiscretizationError,
                      std::map< std::string, walberla::int64_t >&          sqlIntegerProperties,
                      std::map< std::string, double >&                     sqlRealProperties,
                      std::map< std::string, std::string >&                sqlStringProperties,
                      std::map< uint_t, std::map< std::string, double > >& sqlRealPropertiesMG )
{
   WALBERLA_UNUSED( sqlStringProperties );

   StokesFunction u( "u", storage, minLevel, maxLevel );
   StokesFunction f( "f", storage, minLevel, maxLevel );

   StokesFunction uExact( "uExact", storage, minLevel, maxLevel );
   StokesFunction residual( "residual", storage, minLevel, maxLevel );
   StokesFunction error( "error", storage, minLevel, maxLevel );
   StokesFunction tmp( "tmp", storage, minLevel, maxLevel );

   StokesOperator A( storage, minLevel, maxLevel );
   MassOperator   M( storage, minLevel, maxLevel );

   real_t l2ErrorU;
   real_t l2ErrorV;
   real_t l2ErrorP;
   real_t l2ResidualU;
   real_t l2ResidualV;
   real_t l2ResidualP;

   ////////////////////
   // Initialize VTK //
   ////////////////////

   VTKOutput vtkOutput( "vtk", "P2MultigridStokes", storage );
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
      u.u.interpolate( bcU, level, DirichletBoundary );
      u.v.interpolate( exactV, level, DirichletBoundary );

      uExact.u.interpolate( exactU, level, All );
      uExact.v.interpolate( exactV, level, All );
      uExact.p.interpolate( exactP, level, All );
      vertexdof::projectMean( uExact.p, level );

      tmp.u.interpolate( rhsU, level, All );
      tmp.v.interpolate( rhsV, level, All );
      M.apply( tmp.u, f.u, level, All );
      M.apply( tmp.v, f.v, level, All );
   }

   /////////////////////////
   // Misc setup and info //
   /////////////////////////

   real_t discretizationErrorU;
   real_t discretizationErrorV;
   real_t discretizationErrorP;
   if ( calcDiscretizationError )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "l2 discretization error ( u | v | p ) per level:" );
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         calculateDiscretizationErrorStokes< StokesFunction, StokesOperator, MassOperator >(
             storage, level, discretizationErrorU, discretizationErrorV, discretizationErrorP );
         WALBERLA_LOG_INFO_ON_ROOT( "  level " << std::setw( 2 ) << level << ": " << std::scientific << discretizationErrorU
                                               << " | " << discretizationErrorV << " | " << discretizationErrorP );
         sqlRealProperties["l2_discr_error_u_level_" + std::to_string( level )] = discretizationErrorU;
         sqlRealProperties["l2_discr_error_v_level_" + std::to_string( level )] = discretizationErrorV;
         sqlRealProperties["l2_discr_error_p_level_" + std::to_string( level )] = discretizationErrorP;
      }
      WALBERLA_LOG_DEVEL_ON_ROOT( "" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Number of unknowns (including boundary):" )
   uint_t totalDoFs = 0;
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      const uint_t dofsThisLevel = numberOfGlobalDoFs< typename StokesFunction::Tag >( *storage, level );
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
   calculateErrorAndResidualStokes(
       maxLevel, A, u, f, uExact, error, residual, tmp, l2ErrorU, l2ErrorV, l2ErrorP, l2ResidualU, l2ResidualV, l2ResidualP );
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

   real_t avgl2ErrorConvergenceRateU    = 0;
   real_t avgl2ResidualConvergenceRateU = 0;
   real_t avgl2ErrorConvergenceRateP    = 0;
   real_t avgl2ResidualConvergenceRateP = 0;

   real_t l2ErrorReductionU    = 0;
   real_t l2ResidualReductionU = 0;
   real_t l2ErrorReductionP    = 0;
   real_t l2ResidualReductionP = 0;

   sqlRealPropertiesMG[0]["lowercase_l2_error_u"]              = l2ErrorU;
   sqlRealPropertiesMG[0]["lowercase_l2_error_reduction_u"]    = l2ErrorReductionU;
   sqlRealPropertiesMG[0]["lowercase_l2_residual_u"]           = l2ResidualU;
   sqlRealPropertiesMG[0]["lowercase_l2_residual_reduction_u"] = l2ResidualReductionU;

   sqlRealPropertiesMG[0]["lowercase_l2_error_p"]              = l2ErrorP;
   sqlRealPropertiesMG[0]["lowercase_l2_error_reduction_p"]    = l2ErrorReductionP;
   sqlRealPropertiesMG[0]["lowercase_l2_residual_p"]           = l2ResidualP;
   sqlRealPropertiesMG[0]["lowercase_l2_residual_reduction_p"] = l2ResidualReductionP;

   ///////////
   // Solve //
   ///////////

   auto smoother = std::make_shared< UzawaSmoother< StokesOperator > >( storage, minLevel, maxLevel, false, sorRelax );

   const uint_t preconditionerCGIterations = 0;

   auto cgVelocity = std::make_shared< CGSolver< typename StokesOperator::VelocityOperator_T > >(
       storage, minLevel, maxLevel, preconditionerCGIterations, 1e-14 );
   auto preconditioner = std::make_shared< StokesBlockDiagonalPreconditioner< StokesOperator, P1LumpedInvMassOperator > >(
       storage, minLevel, maxLevel, 1, cgVelocity );
   auto coarseGridSolver = std::make_shared< MinResSolver< StokesOperator > >(
       storage, minLevel, maxLevel, coarseGridMaxIterations, 1e-16, preconditioner );

#ifdef HHG_BUILD_WITH_PETSC
   const uint_t petscLevel  = ( numCycles == 0 ? maxLevel : minLevel );
   auto         petscSolver = std::make_shared< PETScLUSolver< StokesOperator > >( storage, petscLevel );
#endif

   auto prolongationOperator = std::make_shared< Prolongation >();
   auto restrictionOperator  = std::make_shared< Restriction >( projectPressureAfterRestriction );

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >( storage,
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
      real_t _l2ErrorU, _l2ErrorV, _l2ErrorP, _l2ResidualU, _l2ResidualV, _l2ResidualP;
      calculateErrorAndResidualStokes( currentLevel,
                                       A,
                                       u,
                                       f,
                                       uExact,
                                       error,
                                       residual,
                                       tmp,
                                       _l2ErrorU,
                                       _l2ErrorV,
                                       _l2ErrorP,
                                       _l2ResidualU,
                                       _l2ResidualV,
                                       _l2ResidualP );
      sqlRealProperties["fmg_l2_error_level_" + std::to_string( currentLevel )] = _l2ErrorU;
      WALBERLA_LOG_INFO_ON_ROOT( "    fmg level " << currentLevel << ": l2 error: " << std::scientific << _l2ErrorU );
   };

   FullMultigridSolver< StokesOperator > fullMultigridSolver(
       storage, multigridSolver, fmgProlongation, minLevel, maxLevel, fmgInnerCycles, postCycle );

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
      calculateErrorAndResidualStokes(
          maxLevel, A, u, f, uExact, error, residual, tmp, l2ErrorU, l2ErrorV, l2ErrorP, l2ResidualU, l2ResidualV, l2ResidualP );
      vtkOutput.write( maxLevel, 1 );
      WALBERLA_LOG_INFO_ON_ROOT( std::setw( 15 ) << 1 << " || " << std::scientific << l2ErrorU << " | " << l2ErrorP << " | "
                                                 << "      " << l2ErrorReductionU << " || " << l2ResidualU << " | " << l2ResidualP
                                                 << " |          " << l2ResidualReductionU << " || " << std::fixed
                                                 << std::setprecision( 2 ) << std::setw( 14 ) << timeCycle << " | "
                                                 << std::setw( 26 ) << timeError << " | " << std::setw( 12 ) << timeVTK << " |" );
   }

   for ( uint_t cycle = 1; cycle <= numCycles; cycle++ )
   {
      const real_t lastl2ErrorU    = l2ErrorU;
      const real_t lastl2ResidualU = l2ResidualU;

      const real_t lastl2ErrorP    = l2ErrorP;
      const real_t lastl2ResidualP = l2ResidualP;

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

      if ( !NEUMANN_PROBLEM )
         vertexdof::projectMean( u.p, maxLevel );

      timer.reset();
      calculateErrorAndResidualStokes(
          maxLevel, A, u, f, uExact, error, residual, tmp, l2ErrorU, l2ErrorV, l2ErrorP, l2ResidualU, l2ResidualV, l2ResidualP );
      timer.end();
      timeError = timer.last();

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

      sqlRealPropertiesMG[cycle]["lowercase_l2_error_u"]              = l2ErrorU;
      sqlRealPropertiesMG[cycle]["lowercase_l2_error_reduction_u"]    = l2ErrorReductionU;
      sqlRealPropertiesMG[cycle]["lowercase_l2_residual_u"]           = l2ResidualU;
      sqlRealPropertiesMG[cycle]["lowercase_l2_residual_reduction_u"] = l2ResidualReductionU;

      sqlRealPropertiesMG[cycle]["lowercase_l2_error_p"]              = l2ErrorP;
      sqlRealPropertiesMG[cycle]["lowercase_l2_error_reduction_p"]    = l2ErrorReductionP;
      sqlRealPropertiesMG[cycle]["lowercase_l2_residual_p"]           = l2ResidualP;
      sqlRealPropertiesMG[cycle]["lowercase_l2_residual_reduction_p"] = l2ResidualReductionP;

      if ( l2ResidualU < L2residualTolerance )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "l2 residual (u) dropped below tolerance." )
         break;
      }
   }

   avgl2ErrorConvergenceRateU /= real_c( numCycles - skipCyclesForAvgConvRate );
   avgl2ResidualConvergenceRateU /= real_c( numCycles - skipCyclesForAvgConvRate );

   avgl2ErrorConvergenceRateP /= real_c( numCycles - skipCyclesForAvgConvRate );
   avgl2ResidualConvergenceRateP /= real_c( numCycles - skipCyclesForAvgConvRate );

   sqlRealProperties["avg_lowercase_l2_error_conv_rate_u"]    = avgl2ErrorConvergenceRateU;
   sqlRealProperties["avg_lowercase_l2_residual_conv_rate_u"] = avgl2ResidualConvergenceRateU;

   sqlRealProperties["avg_lowercase_l2_error_conv_rate_p"]    = avgl2ErrorConvergenceRateP;
   sqlRealProperties["avg_lowercase_l2_residual_conv_rate_p"] = avgl2ResidualConvergenceRateP;

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
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

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
   const uint_t      numProcesses                    = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
   const uint_t      numFacesPerSide                 = mainConf.getParameter< uint_t >( "numFacesPerSide" );
   const std::string discretization                  = mainConf.getParameter< std::string >( "discretization" );
   const uint_t      numCycles                       = mainConf.getParameter< uint_t >( "numCycles" );
   const std::string cycleTypeString                 = mainConf.getParameter< std::string >( "cycleType" );
   const uint_t      fmgInnerCycles                  = mainConf.getParameter< uint_t >( "fmgInnerCycles" );
   const real_t      L2residualTolerance             = mainConf.getParameter< real_t >( "L2residualTolerance" );
   const real_t      sorRelax                        = mainConf.getParameter< real_t >( "sorRelax" );
   const uint_t      preSmoothingSteps               = mainConf.getParameter< uint_t >( "preSmoothingSteps" );
   const uint_t      postSmoothingSteps              = mainConf.getParameter< uint_t >( "postSmoothingSteps" );
   const uint_t      smoothingIncrement              = mainConf.getParameter< uint_t >( "smoothingIncrement" );
   const bool        projectPressureAfterRestriction = mainConf.getParameter< bool >( "projectPressureAfterRestriction" );
   const uint_t      minLevel                        = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel                        = ( discretization == "P1" ? mainConf.getParameter< uint_t >( "maxLevel" ) :
                                                      mainConf.getParameter< uint_t >( "maxLevel" ) - 1 );
   const bool        calculateDiscretizationError    = mainConf.getParameter< bool >( "calculateDiscretizationError" );
   const uint_t      coarseGridMaxIterations         = mainConf.getParameter< uint_t >( "coarseGridMaxIterations" );
   const bool        outputVTK                       = mainConf.getParameter< bool >( "outputVTK" );
   const bool        outputTiming                    = mainConf.getParameter< bool >( "outputTiming" );
   const bool        outputTimingJSON                = mainConf.getParameter< bool >( "outputTimingJSON" );
   const bool        outputSQL                       = mainConf.getParameter< bool >( "outputSQL" );
   const std::string sqlTag                          = mainConf.getParameter< std::string >( "sqlTag", "default" );
   const uint_t      skipCyclesForAvgConvRate        = mainConf.getParameter< uint_t >( "skipCyclesForAvgConvRate" );
   const std::string meshLayout                      = mainConf.getParameter< std::string >( "meshLayout" );

   // parameter checks
   WALBERLA_CHECK( equation == "stokes" || cycleTypeString == "poisson" );
   WALBERLA_CHECK( discretization == "P1" || discretization == "P2" );
   WALBERLA_CHECK( cycleTypeString == "V" || cycleTypeString == "W" );

   const CycleType cycleType = ( cycleTypeString == "V" ? CycleType::VCYCLE : CycleType::WCYCLE );

   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - equation:                      " << equation );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num processes:                 " << numProcesses );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num faces per side:            " << numFacesPerSide );
   WALBERLA_LOG_INFO_ON_ROOT( "  - discretization:                " << discretization );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num cycles:                    " << numCycles );
   WALBERLA_LOG_INFO_ON_ROOT( "  - cycle type:                    " << cycleTypeString );
   WALBERLA_LOG_INFO_ON_ROOT(
       "  - full multigrid:                "
       << ( fmgInnerCycles == 0 ? "no" : "yes, inner cycles per level: " + std::to_string( fmgInnerCycles ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual tolerance:         " << L2residualTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "  - SOR relax:                     " << sorRelax );
   WALBERLA_LOG_INFO_ON_ROOT( "  - pre- / post- / incr-smoothing:         " << preSmoothingSteps << " / " << postSmoothingSteps
                                                                            << " / " << smoothingIncrement );
   WALBERLA_LOG_INFO_ON_ROOT( "  - min / max level:               " << minLevel << " / " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "  - project pressure after restriction: " << ( projectPressureAfterRestriction ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - calculate discretization error: " << ( calculateDiscretizationError ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - coarse grid max itertions (stokes only): " << coarseGridMaxIterations );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output VTK:                    " << ( outputVTK ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing:                 " << ( outputTiming ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing JSON:            " << ( outputTimingJSON ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output SQL:                    " << ( outputSQL ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - SQL tag:                       " << sqlTag );
   WALBERLA_LOG_INFO_ON_ROOT( "  - skip cycles for avg conv rate: " << skipCyclesForAvgConvRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - mesh layout:                   " << meshLayout );
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

   sqlIntegerProperties["num_processes"]              = int64_c( numProcesses );
   sqlIntegerProperties["num_faces_per_side"]         = int64_c( numFacesPerSide );
   sqlStringProperties["discretization"]              = discretization;
   sqlIntegerProperties["num_cycles"]                 = int64_c( numCycles );
   sqlStringProperties["cycle_type"]                  = cycleTypeString;
   sqlIntegerProperties["fmgInnerCycles"]             = int64_c( fmgInnerCycles );
   sqlRealProperties["sor_relax"]                     = sorRelax;
   sqlIntegerProperties["pre_smoothing"]              = int64_c( preSmoothingSteps );
   sqlIntegerProperties["post_smoothing"]             = int64_c( postSmoothingSteps );
   sqlIntegerProperties["incr_smoothing"]             = int64_c( smoothingIncrement );
   sqlIntegerProperties["min_level"]                  = int64_c( minLevel );
   sqlIntegerProperties["max_level"]                  = int64_c( maxLevel );
   sqlIntegerProperties["coarse_grid_max_iterations"] = int64_c( coarseGridMaxIterations );
   sqlIntegerProperties["project_after_restriction"]  = int64_c( projectPressureAfterRestriction );

   ////////////
   // Domain //
   ////////////

   Point2D leftBottom( {0, 0} );
   if ( equation == "stokes" && ( NEUMANN_PROBLEM || COLLIDING_FLOW ) )
   {
      leftBottom = Point2D( {-1, -1} );
   }

   MeshInfo::meshFlavour meshFlavour;
   if ( meshLayout == "CRISS" )
      meshFlavour = MeshInfo::CRISS;
   else if ( meshLayout == "CRISSCROSS" )
      meshFlavour = MeshInfo::CRISSCROSS;
   else
      WALBERLA_ABORT( "Invalid mesh layout." );

   const auto meshInfo = MeshInfo::meshRectangle( leftBottom, Point2D( {1, 1} ), meshFlavour, numFacesPerSide, numFacesPerSide );
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

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   sqlIntegerProperties["num_macro_vertices"] = int64_c( setupStorage.getNumberOfVertices() );
   sqlIntegerProperties["num_macro_edges"]    = int64_c( setupStorage.getNumberOfEdges() );
   sqlIntegerProperties["num_macro_faces"]    = int64_c( setupStorage.getNumberOfFaces() );

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
                                                                preSmoothingSteps,
                                                                postSmoothingSteps,
                                                                smoothingIncrement,
                                                                projectPressureAfterRestriction,
                                                                coarseGridMaxIterations,
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
                                                                preSmoothingSteps,
                                                                postSmoothingSteps,
                                                                smoothingIncrement,
                                                                projectPressureAfterRestriction,
                                                                coarseGridMaxIterations,
                                                                outputVTK,
                                                                skipCyclesForAvgConvRate,
                                                                calculateDiscretizationError,
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

   if ( outputTimingJSON )
   {
      nlohmann::json ttJson;
      walberla::timing::to_json( ttJson, tt );
      std::ofstream jsonOutput;
      jsonOutput.open( "MultigridStudies.json" );
      jsonOutput << ttJson.dump( 4 );
      jsonOutput.close();
   }

   if ( outputSQL )
   {
      WALBERLA_ROOT_SECTION()
      {
         const std::string                  dbFile = "MultigridStudies.db";
         walberla::postprocessing::SQLiteDB db( dbFile );
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
      }
   }
}

} // namespace hhg

int main( int argc, char** argv )
{
   hhg::setup( argc, argv );
}
