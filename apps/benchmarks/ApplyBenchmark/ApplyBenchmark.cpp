
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"
#include "core/math/Constants.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_t;

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;


   auto cfg = std::make_shared< walberla::config::Config >();
   if( env.config() == nullptr )
   {
      auto defaultFile = "./ApplyBenchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   } else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf     = cfg->getBlock( "Parameters" );
   const uint_t                        level        = mainConf.getParameter< uint_t >( "level" );
   std::string                         meshFileName = mainConf.getParameter< std::string >( "mesh" );

   hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );

   //hhg::MeshInfo::meshUnitSquare( 2 );

   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hhg::loadbalancing::roundRobin( setupStorage );

   std::function< real_t( const hhg::Point3D& ) > ones  = []( const hhg::Point3D& ) { return 1.0; };
   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& xx ) {
      //return 5.0;
      return std::sin( walberla::math::M_PI * xx[0] ) + std::cos( walberla::math::M_PI * xx[1] );
      //return ( real_c(std::rand()) / real_c(RAND_MAX));
   };

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage, timingTree );

   if( mainConf.getParameter< bool >( "writeDomain" ) )
   {
      hhg::writeDomainPartitioningVTK( storage, "./output", "domain" );
   }

   hhg::P1Function< double > oneFunc( "x", storage, level, level );
   oneFunc.interpolate( ones, level );
   hhg::P1Function< double > x( "x", storage, level, level );
   hhg::P1Function< double > y( "y", storage, level, level );
   hhg::P1Function< double > z( "z", storage, level, level );
   x.interpolate( exact, level, hhg::Inner );
   //hhg::communication::syncFunctionBetweenPrimitives(x,level);
   hhg::P1ConstantLaplaceOperator mass( storage, level, level );

   const uint_t localDoFs = hhg::numberOfLocalDoFs< hhg::P1FunctionTag >( *storage, level );
   const uint_t totalDoFs = hhg::numberOfGlobalDoFs< hhg::P1FunctionTag >( *storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( "totalDoFs: " << totalDoFs );
   WALBERLA_LOG_INFO( "localDoFs: " << localDoFs << " totalDoFs: " << totalDoFs );

   walberla::WcTimer timer;

   timer.reset();
   LIKWID_MARKER_START( "HyTeG-apply" );
   mass.apply( x, y, level, hhg::Inner );
   LIKWID_MARKER_STOP( "HyTeG-apply" );
   timer.end();
   double hyteg_apply = timer.last();
   WALBERLA_LOG_INFO_ON_ROOT( "HyTeG apply runtime: " << hyteg_apply );

   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "writing VTK output" );
      hhg::VTKOutput vtkOutput("./output", "ApplyBenchmark", storage);
      vtkOutput.add( x );
      vtkOutput.add( z );
      vtkOutput.add( y );
      vtkOutput.write( level );
   }

   walberla::WcTimingTree tt  = timingTree->getReduced();
   auto                   tt2 = tt.getCopyWithRemainder();

   if( mainConf.getParameter< bool >( "printTiming" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( tt2 );
   }

   if( mainConf.getParameter< bool >( "writeJSON" ) )
   {
      nlohmann::json ttjson = nlohmann::json( tt2 );
      std::ofstream  o( "ApplyBenchmarkOutput.json" );
      o << ttjson;
      o.close();
   }

   WALBERLA_LOG_INFO_ON_ROOT( std::scientific <<  " | " << meshFileName << " | " << level << " | " << totalDoFs << " | " << hyteg_apply << " | " );

   LIKWID_MARKER_CLOSE;
}
