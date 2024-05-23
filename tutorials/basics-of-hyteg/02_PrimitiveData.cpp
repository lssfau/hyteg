/*
 * Copyright (c) 2017-2019 Dominik Bartuschat, Dominik Thoennes, Nils Kohl.
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
#include "core/debug/CheckFunctions.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

/**
 * \page 02_PrimitiveData Tutorial 02 - Adding data to primitives
 *
 * \dontinclude tutorials/basics-of-hyteg/02_PrimitiveData.cpp
 *
 * \brief In this tutorial we will add data to primitives
 *
 * \section T02-PrimitiveData-intro Introduction
 *
 * To decouple the simulation from the domain decomposition and communication we implement a
 * generic mechanism to add data structures to primitives.
 *
 * With the approach we can add arbitrary data structures to vertices, edges etc. These could be
 * STL or custom classes.
 *
 * \section T02-PrimitiveData-customclass A custom data structure
 *
 * Typically, we want our primitives to carry simulation data in form of floating point arrays.
 * Therefore let's create a simple class that wraps such a structure.
 *
 * \snippet tutorials/basics-of-hyteg/02_PrimitiveData.cpp SimulationData
 *
 * Since our domain has no idea how to call the constructor, we need to implement a data handling
 * for our custom data structure:
 *
 * \snippet tutorials/basics-of-hyteg/02_PrimitiveData.cpp SimulationDataHandling
 *
 * Our custom data handling must subclass PrimitiveDataHandling, which is templated with
 * the type of the data and the type of the primitive this data belongs to. Since our data
 * shall be the same for all primitives, we type it with Primitive. We could also create special
 * data structures only for vertices or edges. In this case our data handling needs to be typed
 * with the respective subclass of Primitive (e.g. with Vertex or Edge).
 *
 * Since SimulationDataHandling inherits from PrimitiveDataHandling it must implement a few virtual methods.
 * The most important one is PrimitiveDataHandling::initialize(). This function returns a shared pointer to
 * an initialized object. This way, the PrimitiveStorage does not need to know what type of data it stores.
 *
 * \snippet tutorials/basics-of-hyteg/02_PrimitiveData.cpp SimulationDataHandlingInitialize
 *
 * The PrimitiveDataHandling::serialize() and PrimitiveDataHandling::deserialize() methods can be used
 * for runtime load balancing or checkpoint-restart features (but are not required for this tutorial).
 * If you are sure you do not need them, it is possible to inherit from convenience classes like OnlyInitializeDataHandling
 * which implements empty serialize and deserialize methods.
 *
 * \section T02-PrimitiveData-adding Adding the data to the primitives
 *
 * Now we will add and initialize the data. Since it is possible to add different types of data at the same primitives,
 * we need a mechanism to identify our data instance. Therefore we introduce the PrimitiveDataID. It is a templated identifier
 * that allows for type-safe data retrieval.
 *
 * Data is added via the PrimitiveStorage:
 *
 * \snippet tutorials/basics-of-hyteg/02_PrimitiveData.cpp AddingData
 *
 * In this step, each Primitive initializes one instance of SimulationData.
 *
 * \section T02-PrimitiveData-retrieval Data retrieval
 *
 * To obtain the data from a Primitive, simply call its respective method using the PrimitiveDataID:
 *
 * \snippet tutorials/basics-of-hyteg/02_PrimitiveData.cpp DataRetrieval
 *
 * \section T02-PrimitiveData-code Complete Program
 *
 * \include tutorials/basics-of-hyteg/02_PrimitiveData.cpp
 *
 */


/// [SimulationData]
class SimulationData
{
public:

  SimulationData( const uint_t & size ) :
    size_( size ), data_( new real_t[ size ] )
  {}

  // For convenience, let's have everything public
  uint_t   size_;
  real_t * data_;

};
/// [SimulationData]


/// [SimulationDataHandling]
class SimulationDataHandling : public PrimitiveDataHandling< SimulationData, Primitive >
/// [SimulationDataHandling]
{
public:

  ~SimulationDataHandling() {}

  /// [SimulationDataHandlingInitialize]
  std::shared_ptr< SimulationData > initialize( const Primitive * const primitive ) const
  {
    WALBERLA_UNUSED( primitive );
    return std::make_shared< SimulationData >( 42 );
  }
  /// [SimulationDataHandlingInitialize]

  void serialize( const Primitive * const, const PrimitiveDataID< SimulationData, Primitive > &, SendBuffer & ) const
  {
    // write data to buffer
  }

  void deserialize( const Primitive * const, const PrimitiveDataID< SimulationData, Primitive > &, RecvBuffer & ) const
  {
    // read data from buffer
  }
};


class VertexSimulationDataHandling : public OnlyInitializeDataHandling< SimulationData, Vertex >
{
public:

  std::shared_ptr< SimulationData > initialize( const Vertex * const primitive ) const
  {
    WALBERLA_UNUSED( primitive );
    return std::make_shared< SimulationData >( 4711 );
  }

};


void PrimitiveStorageTutorial()
{

  //------------------------------------//
  // From the PrimtivieStorage tutorial //
  //------------------------------------//

  uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  hyteg::MeshInfo meshInfo = MeshInfo::fromGmshFile( "../data/meshes/tri_2el.msh" );
  hyteg::SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );

  // Let's have a debug print
  WALBERLA_LOG_INFO_ON_ROOT( setupStorage );

  hyteg::PrimitiveStorage storage( setupStorage );

  //-----------------//
  // New stuff below //
  //-----------------//

  /// [AddingData]
  // Adding some data to all primitives
  auto simulationDataHandling = std::make_shared< SimulationDataHandling  >();
  PrimitiveDataID< SimulationData, Primitive > simulationDataID;
  storage.addPrimitiveData( simulationDataID, simulationDataHandling, "simulation data" );
  /// [AddingData]

  // Adding some data only to vertices
  auto vertexSimulationDataHandling = std::make_shared< VertexSimulationDataHandling  >();
  PrimitiveDataID< SimulationData, Vertex > vertexSimulationDataID;
  storage.addVertexData( vertexSimulationDataID, vertexSimulationDataHandling, "simulation data (vertices)" );

  /// [DataRetrieval]
  // Gather all the primitives
  PrimitiveStorage::PrimitiveMap primitiveMap;
  storage.getPrimitives( primitiveMap );

  for ( const auto & primitive : primitiveMap )
  {
    WALBERLA_LOG_INFO( "Checking data from Primitive with ID: " << primitive.first );

    // Getting the data via the respective ID
    SimulationData * data = primitive.second->getData( simulationDataID );

    WALBERLA_CHECK_EQUAL( data->size_, 42 );
  }
  /// [DataRetrieval]

  // For nicer output
  WALBERLA_MPI_BARRIER();

  // Check data of the vertices
  for ( const auto & it : storage.getVertices() )
  {
      WALBERLA_LOG_INFO( "Checking data from Vertex with ID: " << it.first );

      // Getting the data via the respective ID
      SimulationData * data = it.second->getData( vertexSimulationDataID );

      WALBERLA_CHECK_EQUAL( data->size_, 4711 );
  }

}

} // namespace hyteg

int main( int argc, char** argv )
{
  walberla::mpi::Environment env( argc, argv );
  walberla::mpi::MPIManager::instance()->useWorldComm();
  hyteg::PrimitiveStorageTutorial();
  return 0;
}

