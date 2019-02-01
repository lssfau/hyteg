
#include "core/Environment.h"
#include "core/config/Config.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
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

void calculateErrorAndResidual( const uint_t&                    level,
                                const P2ConstantLaplaceOperator& A,
                                const P2ConstantMassOperator&    M,
                                const P2Function< real_t >&      u,
                                const P2Function< real_t >&      f,
                                const P2Function< real_t >&      uExact,
                                const P2Function< real_t >&      error,
                                const P2Function< real_t >&      residual,
                                const P2Function< real_t >&      tmp,
                                real_t&                          l2Error,
                                real_t&                          L2Error,
                                real_t&                          l2Residual,
                                real_t&                          L2Residual )
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

void P2MultigridLaplace( const std::shared_ptr< PrimitiveStorage >& storage,
                         const uint_t&                              minLevel,
                         const uint_t&                              maxLevel,
                         const uint_t&                              numVCycles,
                         const uint_t&                              preSmoothingSteps,
                         const uint_t&                              postSmoothingSteps,
                         const bool&                                outputVTK,
                         const uint_t&                              skipCyclesForAvgConvRate )
{
   P2Function< real_t > u( "u", storage, minLevel, maxLevel );
   P2Function< real_t > f( "f", storage, minLevel, maxLevel );

   P2Function< real_t > uExact( "uExact", storage, minLevel, maxLevel );
   P2Function< real_t > residual( "residual", storage, minLevel, maxLevel );
   P2Function< real_t > error( "error", storage, minLevel, maxLevel );
   P2Function< real_t > tmp( "tmp", storage, minLevel, maxLevel );

   P2ConstantLaplaceOperator A( storage, minLevel, maxLevel );
   P2ConstantMassOperator    M( storage, minLevel, maxLevel );

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

   ////////////////
   // Misc setup //
   ////////////////

   calculateErrorAndResidual( maxLevel, A, M, u, f, uExact, error, residual, tmp, l2Error, L2Error, l2Residual, L2Residual );

   WALBERLA_LOG_INFO_ON_ROOT(
       " After vCycle... ||     l2 error |     L2 error | L2 error reduction ||  l2 residual |  L2 residual | L2 residual reduction |" );
   WALBERLA_LOG_INFO_ON_ROOT(
       " ----------------++--------------+--------------+--------------------++--------------+--------------+-----------------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "         initial || " << std::scientific << l2Error << " | " << L2Error << " | "
                                                     << "               --- || " << l2Residual << " | " << L2Residual
                                                     << " |                   --- |" );

   if ( outputVTK )
   {
      vtkOutput.write( maxLevel, 0 );
   }

   real_t avgL2ErrorConvergenceRate    = 0;
   real_t avgL2ResidualConvergenceRate = 0;

   ///////////
   // Solve //
   ///////////

   auto smoother         = std::make_shared< GaussSeidelSmoother< P2ConstantLaplaceOperator > >();
   auto coarseGridSolver = std::make_shared< CGSolver< P2ConstantLaplaceOperator > >( storage, minLevel, minLevel );

   auto prolongationOperator = std::make_shared< P2toP2QuadraticProlongation >();
   auto restrictionOperator  = std::make_shared< P2toP2QuadraticRestriction >();

   GeometricMultigridSolver< P2ConstantLaplaceOperator > multigridSolver( storage,
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

      multigridSolver.solve( A, u, f, maxLevel );
      calculateErrorAndResidual( maxLevel, A, M, u, f, uExact, error, residual, tmp, l2Error, L2Error, l2Residual, L2Residual );

      const real_t L2ErrorReduction    = L2Error / lastL2Error;
      const real_t L2ResidualReduction = L2Residual / lastL2Residual;

      WALBERLA_LOG_INFO_ON_ROOT( std::setw(16) << cycle << " || " << std::scientific << l2Error << " | " << L2Error << " | "
                                                        << "      " << L2ErrorReduction << " || " << l2Residual << " | "
                                                        << L2Residual << " |          " << L2ResidualReduction << " |" );
      if ( outputVTK )
      {
         vtkOutput.write( maxLevel, cycle );
      }

      if ( cycle > skipCyclesForAvgConvRate )
      {
         avgL2ErrorConvergenceRate += L2ErrorReduction;
         avgL2ResidualConvergenceRate += L2ResidualReduction;
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

   const uint_t      numProcesses       = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
   const uint_t      numFacesPerSide    = mainConf.getParameter< uint_t >( "numFacesPerSide" );
   const std::string discretization     = mainConf.getParameter< std::string >( "discretization" );
   const uint_t      numVCycles         = mainConf.getParameter< uint_t >( "numVCycles" );
   const uint_t      preSmoothingSteps  = mainConf.getParameter< uint_t >( "preSmoothingSteps" );
   const uint_t      postSmoothingSteps = mainConf.getParameter< uint_t >( "postSmoothingSteps" );
   const uint_t      minLevel           = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel           = mainConf.getParameter< uint_t >( "maxLevel" );
   const bool        outputVTK          = mainConf.getParameter< bool >( "outputVTK" );
   const uint_t      skipCyclesForAvgConvRate = mainConf.getParameter< uint_t >( "skipCyclesForAvgConvRate" );

   // parameter checks
   WALBERLA_CHECK( discretization == "P1" || discretization == "P2" );

   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num processes:                 " << numProcesses );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num faces per side:            " << numFacesPerSide );
   WALBERLA_LOG_INFO_ON_ROOT( "  - discretization:                " << discretization );
   WALBERLA_LOG_INFO_ON_ROOT( "  - num v-cycles:                  " << numVCycles );
   WALBERLA_LOG_INFO_ON_ROOT( "  - pre- / post-smoothing:         " << preSmoothingSteps << " / " << postSmoothingSteps );
   WALBERLA_LOG_INFO_ON_ROOT( "  - min / max level:               " << minLevel << " / " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "  - output VTK:                    " << outputVTK ? "yes" : "no" );
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

   if ( discretization == "P2" )
   {
      P2MultigridLaplace( storage, minLevel, maxLevel, numVCycles, preSmoothingSteps, postSmoothingSteps, outputVTK, skipCyclesForAvgConvRate );
   }
}

} // namespace hhg

int main( int argc, char** argv )
{
   hhg::setup( argc, argv );
}