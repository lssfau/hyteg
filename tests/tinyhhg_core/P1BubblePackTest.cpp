#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"


using namespace hhg::communication;
using namespace hhg;

class FaceP1BubbleFunctionMemoryDataHandling : public OnlyInitializeDataHandling< FaceP1BubbleFunctionMemory, Face >{
public:
    FaceP1BubbleFunctionMemory * initialize(const Face * const face) const
    {
      //face->memory.push_back(new FaceP1BubbleFunctionMemory());
      WALBERLA_UNUSED(face);
      for (size_t level = 1; level <= 5; ++level)
      {
        P1Bubble::getFaceFunctionMemory(*face, 0)->addlevel(level);
      }
      return P1Bubble::getFaceFunctionMemory(*face, 0);
    }
};

class EdgeP1BubbleFunctionMemoryDataHandling : public OnlyInitializeDataHandling< EdgeP1BubbleFunctionMemory, Edge >{
public:
    EdgeP1BubbleFunctionMemory * initialize(const Edge * const edge) const
    {
      WALBERLA_UNUSED(edge);
      EdgeP1BubbleFunctionMemory *mem = new EdgeP1BubbleFunctionMemory();
      for (size_t level = 1; level <= 5; ++level)
      {
        mem->addlevel(level,2);
      }
      return mem;
    }
};

class VertexP1BubbleFunctionMemoryDataHandling : public OnlyInitializeDataHandling< VertexP1BubbleFunctionMemory, Vertex >{
public:
    VertexP1BubbleFunctionMemory * initialize(const Vertex * const vertex) const
    {
      WALBERLA_UNUSED(vertex);
      VertexP1BubbleFunctionMemory *mem = new VertexP1BubbleFunctionMemory();
      for (size_t level = 1; level <= 5; ++level)
      {
        mem->addlevel(level,vertex->edges.size());
      }
      return mem;
    }
};

int main (int argc, char** argv) {
  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  const uint_t maxLevel = 3;

  hhg::MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, 1  );

  std::shared_ptr<hhg::PrimitiveStorage> storage(new hhg::PrimitiveStorage(uint_c(walberla::MPIManager::instance()->rank()) , setupStorage));

  FaceP1BubbleFunctionMemoryDataHandling faceP1BubbleFunctionMemoryDataHandling;
  EdgeP1BubbleFunctionMemoryDataHandling edgeP1BubbleFunctionMemoryDataHandling;
  VertexP1BubbleFunctionMemoryDataHandling vertexP1BubbleFunctionMemoryDataHandling;
  hhg::PrimitiveDataID<FaceP1BubbleFunctionMemory, Face> faceDataID = storage->addFaceData(faceP1BubbleFunctionMemoryDataHandling,"data");
  hhg::PrimitiveDataID<EdgeP1BubbleFunctionMemory, Edge> edgeDataID = storage->addEdgeData(edgeP1BubbleFunctionMemoryDataHandling,"data");
  hhg::PrimitiveDataID<VertexP1BubbleFunctionMemory, Vertex> vertexDataID = storage->addVertexData(vertexP1BubbleFunctionMemoryDataHandling,"data");

  hhg::PrimitiveStorage::PrimitiveMap primitveMap;
  storage->getPrimitives( primitveMap);

  FaceP1BubbleFunctionMemory *data = storage->beginFaces()->second->getData(faceDataID);
  WALBERLA_UNUSED(data);

  P1BubblePackInfo *packInfo = new P1BubblePackInfo(maxLevel,vertexDataID,edgeDataID,faceDataID,storage);



  std::function<real_t(const hhg::Point3D&)> eight = [](const hhg::Point3D&) { return 8; };
  std::function<real_t(const hhg::Point3D&)> nine = [](const hhg::Point3D&) { return 9; };

  hhg::Face *face0 = storage->beginFaces()->second;
//
  walberla::mpi::SendBuffer sb;
  packInfo->packVertexForEdge(storage->beginVertices()->second,storage->beginEdges()->first,sb);
//  hhg::P1BubbleVertex::packData(level,*edge.v0,0,sb);
//  hhg::P1BubbleVertex::packData(level,*edge.v1,0,sb);
//  walberla::mpi::RecvBuffer rb(sb);
//  unpackVertexData(level,edge,memory_id,rb,*edge.v0);

  hhg::P1BubbleFace::interpolate(*face0,0,eight,maxLevel);

  WALBERLA_UNUSED(face0);

}