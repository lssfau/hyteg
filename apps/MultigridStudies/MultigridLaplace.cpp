
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/timing/TimingJSON.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GaussSeidelSmoother.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"

namespace hhg {

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

std::function< real_t( const hhg::Point3D& ) > rhs = []( const hhg::Point3D& ) { return 0; };
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
}

template < typename Function, typename LaplaceOperator, typename MassOperator, typename Restriction, typename Prolongation >
void MultigridLaplace( const std::shared_ptr< PrimitiveStorage >& storage,
                       const uint_t&                              minLevel,
                       const uint_t&                              maxLevel,
                       const uint_t&                              numVCycles,
                       const real_t&                              L2residualTolerance,
                       const uint_t&                              preSmoothingSteps,
                       const uint_t&                              postSmoothingSteps,
                       const bool&                                outputVTK,
                       const uint_t&                              skipCyclesForAvgConvRate )
{
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

   u.interpolate( exact, maxLevel, DirichletBoundary );
   uExact.interpolate( exact, maxLevel, All );

   tmp.interpolate( rhs, maxLevel, All );
   M.apply( tmp, f, maxLevel, All );

   /////////////////////////
   // Misc setup and info //
   /////////////////////////

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

   walberla::WcTimer timer;
   double            timeError;
   double            timeVTK;
   double            timeCycle;

   timer.reset();
   calculateErrorAndResidual( maxLevel, A, M, u, f, uExact, error, residual, tmp, l2Error, L2Error, l2Residual, L2Residual );
   timer.end();
   timeError = timer.last();

   if ( outputVTK )
   {
      timer.reset();
      vtkOutput.write( maxLevel, 0 );
      timer.end();
      timeVTK = timer.last();
   }

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

   ///////////
   // Solve //
   ///////////

   auto smoother         = std::make_shared< GaussSeidelSmoother< LaplaceOperator > >();
   auto coarseGridSolver = std::make_shared< CGSolver< LaplaceOperator > >( storage, minLevel, minLevel );

   auto prolongationOperator = std::make_shared< Prolongation >();
   auto restrictionOperator  = std::make_shared< Restriction >();

   GeometricMultigridSolver< LaplaceOperator > multigridSolver( storage,
                                                                smoother,
                                                                coarseGridSolver,
                                                                restrictionOperator,
                                                                prolongationOperator,
                                                                minLevel,
                                                                maxLevel,
                                                                preSmoothingSteps,
                                                                postSmoothingSteps,
                                                                0 );

   for ( uint_t cycle = 1; cycle <= numVCycles; cycle++ )
   {
      const real_t lastL2Error    = L2Error;
      const real_t lastL2Residual = L2Residual;

      timer.reset();
      multigridSolver.solve( A, u, f, maxLevel );
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

      const real_t L2ErrorReduction    = L2Error / lastL2Error;
      const real_t L2ResidualReduction = L2Residual / lastL2Residual;

      WALBERLA_LOG_INFO_ON_ROOT( std::setw( 15 ) << cycle << " || " << std::scientific << l2Error << " | " << L2Error << " | "
                                                 << "      " << L2ErrorReduction << " || " << l2Residual << " | " << L2Residual
                                                 << " |          " << L2ResidualReduction << " || " << std::fixed
                                                 << std::setprecision( 2 ) << std::setw( 14 ) << timeCycle << " | "
                                                 << std::setw( 26 ) << timeError << " | " << std::setw( 12 ) << timeVTK << " |" );

      if ( cycle > skipCyclesForAvgConvRate )
      {
         avgL2ErrorConvergenceRate += L2ErrorReduction;
         avgL2ResidualConvergenceRate += L2ResidualReduction;
      }

      if ( L2Residual < L2residualTolerance )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "L2 residual dropped below tolerance." )
         break;
      }
   }

   avgL2ErrorConvergenceRate /= real_c( numVCycles - skipCyclesForAvgConvRate );
   avgL2ResidualConvergenceRate /= real_c( numVCycles - skipCyclesForAvgConvRate );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Average convergence rates:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 error:    " << std::scientific << avgL2ErrorConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual: " << std::scientific << avgL2ResidualConvergenceRate );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

void setup( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "///////////////////////" );
   WALBERLA_LOG_INFO_ON_ROOT( "// Multigrid Studies //" );
   WALBERLA_LOG_INFO_ON_ROOT( "///////////////////////" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./MultigridLaplace.prm";
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

   const uint_t      numProcesses             = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
   const uint_t      numFacesPerSide          = mainConf.getParameter< uint_t >( "numFacesPerSide" );
   const std::string discretization           = mainConf.getParameter< std::string >( "discretization" );
   const uint_t      numVCycles               = mainConf.getParameter< uint_t >( "numVCycles" );
   const real_t      L2residualTolerance      = mainConf.getParameter< real_t >( "L2residualTolerance" );
   const uint_t      preSmoothingSteps        = mainConf.getParameter< uint_t >( "preSmoothingSteps" );
   const uint_t      postSmoothingSteps       = mainConf.getParameter< uint_t >( "postSmoothingSteps" );
   const uint_t      minLevel                 = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel                 = mainConf.getParameter< uint_t >( "maxLevel" );
   const bool        outputVTK                = mainConf.getParameter< bool >( "outputVTK" );
   const bool        outputTiming             = mainConf.getParameter< bool >( "outputTiming" );
   const bool        outputTimingJSON         = mainConf.getParameter< bool >( "outputTimingJSON" );
   const uint_t      skipCyclesForAvgConvRate = mainConf.getParameter< uint_t >( "skipCyclesForAvgConvRate" );

   // parameter checks
   WALBERLA_CHECK( discretization == "P1" || discretization == "P2" );

   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num processes:                 " << numProcesses );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num faces per side:            " << numFacesPerSide );
   WALBERLA_LOG_INFO_ON_ROOT( "  - discretization:                " << discretization );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num v-cycles:                  " << numVCycles );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual tolerance:         " << L2residualTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "  - pre- / post-smoothing:         " << preSmoothingSteps << " / " << postSmoothingSteps );
   WALBERLA_LOG_INFO_ON_ROOT( "  - min / max level:               " << minLevel << " / " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output VTK:                    " << ( outputVTK ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing:                 " << ( outputTiming ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing JSON:            " << ( outputTimingJSON ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - skip cycles for avg conv rate: " << skipCyclesForAvgConvRate );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   ////////////
   // Domain //
   ////////////

   const auto meshInfo =
       MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1, 1} ), MeshInfo::CRISS, numFacesPerSide + 1, numFacesPerSide + 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   if ( outputVTK )
   {
      writeDomainPartitioningVTK( storage, "vtk", "Domain" );
   }

   auto globalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( globalInfo );

   if ( discretization == "P1" )
   {
      MultigridLaplace< P1Function< real_t >,
                        P1ConstantLaplaceOperator,
                        P1MassOperator,
                        P1toP1LinearRestriction,
                        P1toP1LinearProlongation >( storage,
                                                    minLevel,
                                                    maxLevel,
                                                    numVCycles,
                                                    L2residualTolerance,
                                                    preSmoothingSteps,
                                                    postSmoothingSteps,
                                                    outputVTK,
                                                    skipCyclesForAvgConvRate );
   }
   else if ( discretization == "P2" )
   {
      MultigridLaplace< P2Function< real_t >,
                        P2ConstantLaplaceOperator,
                        P2ConstantMassOperator,
                        P2toP2QuadraticRestriction,
                        P2toP2QuadraticProlongation >( storage,
                                                       minLevel,
                                                       maxLevel,
                                                       numVCycles,
                                                       L2residualTolerance,
                                                       preSmoothingSteps,
                                                       postSmoothingSteps,
                                                       outputVTK,
                                                       skipCyclesForAvgConvRate );
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
}

} // namespace hhg

int main( int argc, char** argv )
{
   hhg::setup( argc, argv );
}