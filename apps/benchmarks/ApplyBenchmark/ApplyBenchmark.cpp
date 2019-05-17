#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   LIKWID_MARKER_THREADINIT;

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./ApplyBenchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf        = cfg->getBlock( "Parameters" );
   const uint_t                        level           = mainConf.getParameter< uint_t >( "level" );
   const uint_t                        facesPerProcess = mainConf.getParameter< uint_t >( "facesPerProcess" );
   const uint_t                        numProc         = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   hhg::MeshInfo meshInfo = hhg::MeshInfo::meshFaceChain( numProc * facesPerProcess );

   hhg::SetupPrimitiveStorage setupStorage( meshInfo, numProc );
   hhg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hhg::PrimitiveStorage >  storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage, timingTree );

   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& xx ) {
      return std::sin( walberla::math::M_PI * xx[0] ) + std::cos( walberla::math::M_PI * xx[1] );
   };

   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );

   ///// Functions / operators / allocation /////

   hhg::P1Function< double > src( "src", storage, level, level );
   hhg::P1Function< double > dst( "dst", storage, level, level );

   hhg::P1ConstantLaplaceOperator laplace( storage, level, level );

   src.interpolate( exact, level, hhg::Inner );

   const uint_t localDoFs = hhg::numberOfLocalDoFs< hhg::P1FunctionTag >( *storage, level );
   const uint_t totalDoFs = hhg::numberOfGlobalDoFs< hhg::P1FunctionTag >( *storage, level );

   WALBERLA_LOG_INFO( "localDoFs: " << localDoFs << " totalDoFs: " << totalDoFs )

   walberla::WcTimer timer;
   uint_t            iterations = 1;
   do
   {
      timer.reset();
      LIKWID_MARKER_START( "HyTeG-apply" );
      for ( uint_t i = 0; i < iterations; ++i )
      {
         laplace.apply( src, dst, level, hhg::Inner );
      }
      LIKWID_MARKER_STOP( "HyTeG-apply" );
      timer.end();
      iterations *= 2;
   } while ( timer.last() < 1 );

   iterations /= 2;

   double hyteg_apply = timer.last();
   WALBERLA_LOG_INFO_ON_ROOT( "HyTeG apply runtime: " << hyteg_apply )

   walberla::WcTimingTree tt  = timingTree->getReduced();
   auto                   tt2 = tt.getCopyWithRemainder();

   if ( mainConf.getParameter< bool >( "printTiming" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( tt2 )
   }

   if ( mainConf.getParameter< bool >( "writeJSON" ) )
   {
      nlohmann::json ttjson = nlohmann::json( tt2 );
      std::ofstream  o( "ApplyBenchmarkOutput.json" );
      o << ttjson;
      o.close();
   }

   const uint_t globalInnerDoFs = hhg::numberOfGlobalInnerDoFs< hhg::P1FunctionTag >( *storage, level );
   const real_t glups           = real_c( globalInnerDoFs * iterations ) / 1e9 / hyteg_apply;
   const real_t gflops          = real_c( globalInnerDoFs * iterations * 13 ) / 1e9 / hyteg_apply;

   WALBERLA_LOG_INFO_ON_ROOT( hhg::format(
       "%10s|%10s|%10s|%10s|%5s|%5s|%7s", "Time (s)", "GDoF/s", "GFLOP/s", " DoFs ", "Level", "Procs", "face/proc" ) )
   WALBERLA_LOG_INFO_ON_ROOT( hhg::format( "%10.3e|%10.3e|%10.3e|%10.3e|%5u|%5u|%7u",
                                           hyteg_apply,
                                           glups,
                                           gflops,
                                           real_c( globalInnerDoFs ),
                                           level,
                                           numProc,
                                           facesPerProcess ) )
   //WALBERLA_LOG_INFO_ON_ROOT( std::scientific << " | " << level << " | " << totalDoFs << " | " << hyteg_apply << " | " )

   LIKWID_MARKER_CLOSE;
}
