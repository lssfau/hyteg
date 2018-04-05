#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

//using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

//using namespace hhg;

int main( int argc, char** argv )
{
   typedef double real_t;
   /// [Create Environment]
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   /// [Create Environment]

   /// [Get Parameters]
   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "../data/param/GMG_P1.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t      minLevel         = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel         = parameters.getParameter< uint_t >( "maxLevel" );
   const uint_t      max_outer_iter   = parameters.getParameter< uint_t >( "max_outer_iter" );
   const uint_t      max_coarse_iter  = parameters.getParameter< uint_t >( "max_coarse_iter" );
   const real_t      mg_tolerance     = parameters.getParameter< real_t >( "mg_tolerance" );
   const real_t      coarse_tolerance = parameters.getParameter< real_t >( "coarse_tolerance" );
   const std::string meshFile         = parameters.getParameter< std::string >( "mesh" );
   /// [Get Parameters]

   /// [Create Primitive Storage]
   std::shared_ptr< hhg::PrimitiveStorage > storage = hhg::PrimitiveStorage::createFromGmshFile(meshFile);
   /// [Create Primitive Storage]

   /// [Create Functions]
   P1Function< real_t > residuum( "residuum", storage, minLevel, maxLevel );
   /// [Create Functions]
}