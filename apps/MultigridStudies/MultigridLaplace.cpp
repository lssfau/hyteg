
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/timing/TimingJSON.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1QuadraticProlongation.hpp"
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
#include "tinyhhg_core/solvers/FullMultigridSolver.hpp"
#include "tinyhhg_core/solvers/GaussSeidelSmoother.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"

#include "postprocessing/sqlite/SQLite.h"

namespace hhg {

using walberla::int64_c;

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
                       const uint_t&                                          fmgInnerCycles,
                       const real_t&                                        L2residualTolerance,
                       const uint_t&                                        preSmoothingSteps,
                       const uint_t&                                        postSmoothingSteps,
                       const bool&                                          outputVTK,
                       const uint_t&                                        skipCyclesForAvgConvRate,
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

   auto smoother         = std::make_shared< GaussSeidelSmoother< LaplaceOperator > >();
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

   FullMultigridSolver< LaplaceOperator > fullMultigridSolver( storage, multigridSolver, fmgProlongation, minLevel, maxLevel, fmgInnerCycles );

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
   const uint_t      numCycles                = mainConf.getParameter< uint_t >( "numCycles" );
   const std::string cycleTypeString          = mainConf.getParameter< std::string >( "cycleType" );
   const uint_t      fmgInnerCycles           = mainConf.getParameter< uint_t >( "fmgInnerCycles" );
   const real_t      L2residualTolerance      = mainConf.getParameter< real_t >( "L2residualTolerance" );
   const uint_t      preSmoothingSteps        = mainConf.getParameter< uint_t >( "preSmoothingSteps" );
   const uint_t      postSmoothingSteps       = mainConf.getParameter< uint_t >( "postSmoothingSteps" );
   const uint_t      minLevel                 = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel                 = mainConf.getParameter< uint_t >( "maxLevel" );
   const bool        outputVTK                = mainConf.getParameter< bool >( "outputVTK" );
   const bool        outputTiming             = mainConf.getParameter< bool >( "outputTiming" );
   const bool        outputTimingJSON         = mainConf.getParameter< bool >( "outputTimingJSON" );
   const bool        outputSQL                = mainConf.getParameter< bool >( "outputSQL" );
   const std::string sqlTag                   = mainConf.getParameter< std::string >( "sqlTag", "default" );
   const uint_t      skipCyclesForAvgConvRate = mainConf.getParameter< uint_t >( "skipCyclesForAvgConvRate" );

   // parameter checks
   WALBERLA_CHECK( discretization == "P1" || discretization == "P2" );
   WALBERLA_CHECK( cycleTypeString == "V" || cycleTypeString == "W" );

   const CycleType cycleType = ( cycleTypeString == "V" ? CycleType::VCYCLE : CycleType::WCYCLE );

   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num processes:                 " << numProcesses );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num faces per side:            " << numFacesPerSide );
   WALBERLA_LOG_INFO_ON_ROOT( "  - discretization:                " << discretization );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num cycles:                    " << numCycles );
   WALBERLA_LOG_INFO_ON_ROOT( "  - cycle type:                    " << cycleTypeString );
   WALBERLA_LOG_INFO_ON_ROOT( "  - full multigrid:                " << ( fmgInnerCycles == 0 ? "no" : "yes, inner cycles per level: " + std::to_string( fmgInnerCycles ) ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - L2 residual tolerance:         " << L2residualTolerance );
   WALBERLA_LOG_INFO_ON_ROOT( "  - pre- / post-smoothing:         " << preSmoothingSteps << " / " << postSmoothingSteps );
   WALBERLA_LOG_INFO_ON_ROOT( "  - min / max level:               " << minLevel << " / " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output VTK:                    " << ( outputVTK ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing:                 " << ( outputTiming ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output timing JSON:            " << ( outputTimingJSON ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output SQL:                    " << ( outputSQL ? "yes" : "no" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "  - SQL tag:                       " << sqlTag );
   WALBERLA_LOG_INFO_ON_ROOT( "  - skip cycles for avg conv rate: " << skipCyclesForAvgConvRate );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   /////////
   // SQL //
   /////////

   std::map< std::string, walberla::int64_t >          sqlIntegerProperties;
   std::map< std::string, double >                     sqlRealProperties;
   std::map< std::string, std::string >                sqlStringProperties;
   std::map< uint_t, std::map< std::string, double > > sqlRealPropertiesMG;

   sqlStringProperties["tag"] = sqlTag;

   sqlIntegerProperties["num_processes"]      = int64_c( numProcesses );
   sqlIntegerProperties["num_faces_per_side"] = int64_c( numFacesPerSide );
   sqlStringProperties["discretization"]      = discretization;
   sqlIntegerProperties["num_cycles"]         = int64_c( numCycles );
   sqlStringProperties["cycle_type"]          = cycleTypeString;
   sqlIntegerProperties["fmgInnerCycles"]     = int64_c( fmgInnerCycles );
   sqlIntegerProperties["pre_smoothing"]      = int64_c( preSmoothingSteps );
   sqlIntegerProperties["post_smoothing"]     = int64_c( postSmoothingSteps );
   sqlIntegerProperties["min_level"]          = int64_c( minLevel );
   sqlIntegerProperties["max_level"]          = int64_c( maxLevel );

   ////////////
   // Domain //
   ////////////

   const auto meshInfo =
       MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1, 1} ), MeshInfo::CRISS, numFacesPerSide + 1, numFacesPerSide + 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
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

   if ( discretization == "P1" )
   {
      MultigridLaplace< P1Function< real_t >,
                        P1ConstantLaplaceOperator,
                        P1MassOperator,
                        P1toP1LinearRestriction,
                        P1toP1LinearProlongation,
                        P1toP1QuadraticProlongation >( storage,
                                                       minLevel,
                                                       maxLevel,
                                                       numCycles,
                                                       cycleType,
                                                       fmgInnerCycles,
                                                       L2residualTolerance,
                                                       preSmoothingSteps,
                                                       postSmoothingSteps,
                                                       outputVTK,
                                                       skipCyclesForAvgConvRate,
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
                                                       preSmoothingSteps,
                                                       postSmoothingSteps,
                                                       outputVTK,
                                                       skipCyclesForAvgConvRate,
                                                       sqlIntegerProperties,
                                                       sqlRealProperties,
                                                       sqlStringProperties,
                                                       sqlRealPropertiesMG );
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
