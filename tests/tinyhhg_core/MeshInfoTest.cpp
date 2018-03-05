
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/tinyhhg.hpp"

#include "core/Filesystem.h"
#include "boost/algorithm/string/predicate.hpp"

namespace hhg {

static void testMeshInfo()
{
  const std::string meshFileDir = "../../data/meshes/";

  std::vector< std::string > gmshFiles;

  walberla::filesystem::recursive_directory_iterator dirIterator( meshFileDir );
  walberla::filesystem::recursive_directory_iterator dirIteratorEnd;

  for ( ; dirIterator != dirIteratorEnd; dirIterator++ )
  {
    std::string dirOrFilePath = dirIterator->path().string();
    if ( boost::algorithm::ends_with( dirOrFilePath, ".msh" ) )
    {
      gmshFiles.push_back( dirOrFilePath );
    }
  }

  for ( const auto & gmshFile : gmshFiles )
  {
    MeshInfo meshInfo = MeshInfo::fromGmshFile( gmshFile );
    SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
    loadbalancing::roundRobin( setupStorage );
    WALBERLA_LOG_INFO_ON_ROOT( setupStorage );
    PrimitiveStorage storage( setupStorage );
  }
}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testMeshInfo();

   return EXIT_SUCCESS;
}
