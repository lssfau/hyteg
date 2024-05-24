/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "core/mpi/Environment.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hyteg {

/**
 * \page BA.01_PrimitiveStorage Tutorial BA.01 - Setting up a PrimitiveStorage
 *
 * \dontinclude tutorials/basics-for-apps/BA.01_PrimitiveStorage.cpp
 *
 * \brief In this tutorial we will set up a distributed simulation domain using the PrimitiveStorage class.
 *
 * \section T01-intro Introduction
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
 * \section T01-mesh Creating a MeshInfo
 *
 * The MeshInfo class is used to load / parse mesh files of different formats and store the
 * information in an instance. This extra step was introduced to provide an interface for
 * different mesh file formats.
 *
 * To parse a mesh file simple use the respective static method of the class. To parse a
 * file in the GMSH format, run:
 * \snippet{trimleft} this MeshInfo
 *
 * Alternatively you can use one of the inline mesh generators available in HyTeG. The
 * code below e.g. will generate a regular criss-cross mesh composed of (3 x 2) sub-cells
 * (each split into four triangles) on the 2D rectangle [-2,1] x [0,3]:
 *
 * \snippet tutorials/basics-for-apps/99_extraExamples.cpp MeshRectangle
 *
 * Inline generators are provided for a selection of standard geometries such as e.g.
 * rectangle, annulus, thick spherical shell. See the MeshInfo class documentation for
 * a complete list and their options.
 *
 * \section T01-setupstorage Creating a SetupPrimitiveStorage
 *
 * Before we create the distributed storage, we introduce another intermediate step:
 * we create a SetupPrimitiveStorage.
 *
 * This class creates the respective Primitive instances 
 * from the MeshInfo. The resulting storage object is not distributed. However, in this phase the
 * primitives only carry metadata and no actual simulation data is allocated on them.
 *
 * \snippet{trimleft} this SetupPrimitiveStorage
 *
 * The main purpose of the SetupPrimitiveStorage is to assign the primitives to the processes.
 * We achieve this by calling a load balancing function on it. There are some simple load balancers
 * available in the library. However, it will call a default loadbalancer from the library so that
 * the user does not have to call it in the application. Anyway, here is an example of how to use
 * them in practice:
 *
 * \snippet{trimleft} this Loadbalancing
 *
 * \section T01-storage Creating a distributed PrimitiveStorage
 *
 * Finally we are able to create the distributed PrimitiveStorage instance. It requires us to specify
{Paraview/link_002_tomo_SEMUCB_T_angle_}{0}{35}}; * the rank it is located on.
 *
 * \snippet{trimleft} this PrimitiveStorage
 *
 * You can see that the data structure is now distributed since different ranks carry different primitives.
 * Querying the number of primitives should therefore (in general) result in different numbers on different processes:
 *
 * \snippet{trimleft} this PrimitiveStorageNumberOfVertices
 *
 * \section code Complete Program
 * \include tutorials/basics-for-apps/BA.01_PrimitiveStorage.cpp
 *
 *
 *
 */

void PrimitiveStorageTutorial()
{
  uint_t rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );
  uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  /// [MeshInfo]
  hyteg::MeshInfo meshInfo = MeshInfo::fromGmshFile( "../data/meshes/tri_2el.msh" );
  /// [MeshInfo]

  /// [SetupPrimitiveStorage]
  hyteg::SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
  /// [SetupPrimitiveStorage]

  /// [Loadbalancing]
  hyteg::loadbalancing::roundRobin( setupStorage );
  /// [Loadbalancing]

  // Let's have a debug print
  WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

  /// [PrimitiveStorage]
  hyteg::PrimitiveStorage storage( setupStorage );
  /// [PrimitiveStorage]

  // For nicer output
  WALBERLA_MPI_BARRIER();

  /// [PrimitiveStorageNumberOfVertices]
  WALBERLA_LOG_INFO( "Number of vertices on rank " << rank << ": " << storage.getNumberOfLocalVertices() );
  /// [PrimitiveStorageNumberOfVertices]

}

} // namespace hyteg

int main( int argc, char** argv )
{
  walberla::mpi::Environment env( argc, argv );
  walberla::mpi::MPIManager::instance()->useWorldComm();
  hyteg::PrimitiveStorageTutorial();
  return 0;
}
