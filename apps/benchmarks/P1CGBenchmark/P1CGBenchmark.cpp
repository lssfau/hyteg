#include <tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp>
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"
#include "core/timing/TimingNode.h"
#include "core/extern/json.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"

using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   //walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared<walberla::config::Config>();
   if( env.config() == nullptr ) {
      auto defaultFile = "./P1CGBenchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT("No Parameter file given loading default parameter file: " << defaultFile);
      cfg->readParameterFile( defaultFile );
   } else {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t level = mainConf.getParameter< uint_t >( "level" );

   const double      tolerance = 1e-15;
   const uint_t      maxIter   = 1000;

   MeshInfo                            meshInfo = MeshInfo::meshUnitSquare( 2 );
   SetupPrimitiveStorage               setupStorage( meshInfo, 1 );
   //auto storage = PrimitiveStorage( setupStorage );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   //auto storage = PrimitiveStorage::createFromGmshFile( meshFile );
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   storage->setTimingTree(timingTree);

   hhg::P1Function< double > r( "r", storage, level, level );
   hhg::P1Function< double > f( "f", storage, level, level );
   hhg::P1Function< double > u( "u", storage, level, level );
   hhg::P1Function< double > u_exact( "u_exact", storage, level, level );
   hhg::P1Function< double > err( "err", storage, level, level );
   hhg::P1Function< double > npoints_helper( "npoints_helper", storage, level, level );

   hhg::P1MassOperator    M( storage, level, level );
   hhg::P1ConstantLaplaceOperator L( storage, level, level );

   std::function< double( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< double( const hhg::Point3D& ) > rhs = []( const hhg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< double( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hhg::DirichletBoundary );
   u_exact.interpolate( exact, level );
   npoints_helper.interpolate( rhs, level );
   M.apply( npoints_helper, f, level, hhg::All );

   auto solver = hhg::CGSolver< hhg::P1Function< double >, hhg::P1ConstantLaplaceOperator >( storage, level, level );

   solver.solve( L, u, f, r, level, tolerance, maxIter, hhg::Inner, false );

   err.assign( {1.0, -1.0}, {&u, &u_exact}, level );
   npoints_helper.interpolate( ones, level );

   const double npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   const double discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   //WALBERLA_CHECK_LESS( discr_l2_err, 1.2e-5 );


   walberla::WcTimingTree tt = timingTree->getReduced();
   auto tt2 = tt.getCopyWithRemainder();
   nlohmann::json ttjson = nlohmann::json(tt2);
   std::ofstream o("P1CGBenchmarkOutput.json");
   o << ttjson;
   o.close();

//   WALBERLA_LOG_INFO_ON_ROOT( tt );
//   WALBERLA_LOG_INFO_ON_ROOT( tt2 );
//   std::cout << ttjson.dump(2) << std::endl;


   return 0;
}
