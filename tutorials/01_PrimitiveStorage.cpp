#include "core/mpi/Environment.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hhg {

/**
 * \page 01_PrimitiveStorage Setting up a PrimitiveStorage
 *
 * \dontinclude tutorials/01_PrimitiveStorage.cpp
 *
 * \brief In this tutorial we will set up a distributed simulation domain using the PrimitiveStorage class.
 *
 * \section intro Introduction
 *
 * A PrimitiveStorage contains references to all local Primitive instances
 * and therefore implements a totally distributed data structure that represents the simulation domain.
 *
 * It stores information about the neighborhood of the local primitives. Additionally it
 * is used to add simulation data to the primitives and provides iterators to loop over the local primitives.
 *
 * At the moment, a PrimitiveStorage is constructed via a three step approach:
 *
 * -# creating a MeshInfo object from a mesh file
 * -# creating a SetupPrimitiveStorage from the MeshInfo object
 * -# creating a distributed PrimitiveStorage from the SetupPrimitiveStorage
 *
 * \section mesh Creating a MeshInfo
 *
 * The MeshInfo class is used to load / parse mesh files of different formats and store the
 * information in an instance. This extra step was introduced to provide an interface for
 * different mesh file formats.
 *
 * To parse a mesh file simple use the respective static method of the class. To parse a
 * file in the GMSH format, run:
 *
 * \snippet tutorials/01_PrimitiveStorage.cpp MeshInfo
 *
 * \section setupstorage Creating a SetupPrimitiveStorage
 *
 * Before we create the distributed storage, we introduce another intermediate step:
 * we create a SetupPrimitiveStorage.
 *
 * This class creates the respective Primitive instances 
 * from the MeshInfo. The resulting storage object is not distributed. However, in this phase the
 * primitives only carry metadata and no actual simulation data is allocated on them.
 *
 * \snippet tutorials/01_PrimitiveStorage.cpp SetupPrimitiveStorage
 *
 * The main purpose of the SetupPrimitiveStorage is to assign the primitives to the processes.
 * We achieve this by calling a load balancing function on it. There are some simple load balancers
 * available in the library.
 *
 * \snippet tutorials/01_PrimitiveStorage.cpp Loadbalancing
 *
 * \section storage Creating a distributed PrimitiveStorage
 *
 * Finally we are able to create the distributed PrimitiveStorage instance. It requires us to specify
 * the rank it is located on.
 *
 * \snippet tutorials/01_PrimitiveStorage.cpp PrimitiveStorage
 *
 * You can see that the data structure is now distributed since different ranks carry different primitives.
 * Querying the number of primitives should therefore (in general) result in different numbers on different processes:
 *
 * \snippet tutorials/01_PrimitiveStorage.cpp PrimitiveStorageNumberOfVertices
 *
 * \section code Complete Program
 * \include tutorials/01_PrimitiveStorage.cpp
 *
 *
 *
 */

void PrimitiveStorageTutorial()
{
  uint_t rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );
  uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  /// [MeshInfo]
  hhg::MeshInfo meshInfo = MeshInfo::fromGmshFile( "../data/meshes/tri_2el.msh" );
  /// [MeshInfo]

  /// [SetupPrimitiveStorage]
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
  /// [SetupPrimitiveStorage]

  /// [Loadbalancing]
  hhg::loadbalancing::roundRobin( setupStorage );
  /// [Loadbalancing]

  // Let's have a debug print
  WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

  /// [PrimitiveStorage]
  hhg::PrimitiveStorage storage( setupStorage );
  /// [PrimitiveStorage]

  // For nicer output
  WALBERLA_MPI_BARRIER();

  /// [PrimitiveStorageNumberOfVertices]
  WALBERLA_LOG_INFO( "Number of vertices on rank " << rank << ": " << storage.getNumberOfLocalVertices() );
  /// [PrimitiveStorageNumberOfVertices]

}

} // namespace hhg

int main( int argc, char** argv )
{
  walberla::mpi::Environment env( argc, argv );
  walberla::mpi::MPIManager::instance()->useWorldComm();
  hhg::PrimitiveStorageTutorial();
  return 0;
}


