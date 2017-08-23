
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "tinyhhg_core/tinyhhg.hpp"

namespace hhg {

struct TestData
{
  uint_t ownID;
  std::vector< uint_t > neighborIDs;
};

struct TestDataHandling : OnlyInitializeDataHandling< TestData, Primitive >
{
  virtual std::shared_ptr< TestData > initialize( const Primitive * const primitive ) const
  {
    auto data = std::make_shared< TestData >();
    data->ownID = primitive->getID().getID();
    return data;
  }
};

class TestPackInfo : public communication::PackInfo
{
public:

  TestPackInfo( PrimitiveDataID< TestData, Primitive > & dataID ) : dataID_( dataID ) {}

  virtual void packVertexForEdge( const Vertex *sender, const PrimitiveID & /* receiver */, walberla::mpi::SendBuffer & buffer )
  {
    TestData * data = sender->getData( dataID_ );
    buffer << data->ownID;
  }

  virtual void unpackEdgeFromVertex( Edge *receiver, const PrimitiveID & /* sender */, walberla::mpi::RecvBuffer & buffer )
  {
    TestData * data = receiver->getData( dataID_ );
    uint_t vertexData;
    buffer >> vertexData;
    data->neighborIDs.push_back( vertexData );
  }

  virtual void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver)
  {
    TestData * vertexData = sender->getData( dataID_ );
    TestData * edgeData   = receiver->getData( dataID_ );
    edgeData->neighborIDs.push_back( vertexData->ownID );
  }

  // left other methods empty for this tutorial

  virtual void packEdgeForVertex( const Edge *, const PrimitiveID &, walberla::mpi::SendBuffer & ) {}

  virtual void unpackVertexFromEdge( Vertex *, const PrimitiveID &, walberla::mpi::RecvBuffer & ) {}

  virtual void communicateLocalEdgeToVertex( const Edge *, Vertex * ) {}

  virtual void packEdgeForFace( const Edge *, const PrimitiveID &, walberla::mpi::SendBuffer & ) {}

  virtual void unpackFaceFromEdge( Face *, const PrimitiveID &, walberla::mpi::RecvBuffer & ) {}

  virtual void communicateLocalEdgeToFace( const Edge *, Face * ) {}

  virtual void packFaceForEdge( const Face *, const PrimitiveID &, walberla::mpi::SendBuffer & ) {}

  virtual void unpackEdgeFromFace( Edge *, const PrimitiveID &, walberla::mpi::RecvBuffer & ) {}

  virtual void communicateLocalFaceToEdge( const Face *, Edge * ) {}


private:

  PrimitiveDataID< TestData, Primitive > dataID_;

};

}

int main()
{
  uint_t numProcesses = walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( "../data/meshes/tri_2el.msh" );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );

  hhg::loadbalancing::roundRobin( setupStorage );

  auto storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

  hhg::PrimitiveDataID< hhg::TestData, hhg::Primitive > dataID;
  auto testDataHandling = std::make_shared< hhg::TestDataHandling >();
  storage->addPrimitiveData( dataID, testDataHandling, "test data" );

  std::shared_ptr< hhg::TestPackInfo > packInfo = std::make_shared< hhg::TestPackInfo >( dataID );
}



