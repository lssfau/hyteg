#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"


using namespace hhg::communication;
using namespace hhg;

class FaceP1BubbleFunctionMemoryDataHandling : public OnlyInitializeDataHandling< FaceP1BubbleFunctionMemory, Face >{
public:
    FaceP1BubbleFunctionMemory * initialize(const Face * const face) const
    {
      WALBERLA_UNUSED(face);
      FaceP1BubbleFunctionMemory *mem = new FaceP1BubbleFunctionMemory();
      for (size_t level = 1; level <= 5; ++level)
      {
        mem->addlevel(level);
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

  hhg::PrimitiveStorage storage(walberla::MPIManager::instance()->rank() , setupStorage);

  FaceP1BubbleFunctionMemoryDataHandling faceP1BubbleFunctionMemoryDataHandling;
  hhg::PrimitiveDataID<FaceP1BubbleFunctionMemory, Face> dataID = storage.addFaceData(faceP1BubbleFunctionMemoryDataHandling,"data");

  hhg::PrimitiveStorage::PrimitiveMap primitveMap;
  storage.getPrimitives( primitveMap);

  FaceP1BubbleFunctionMemory *data = storage.beginFaces()->second->getData(dataID);


  std::function<real_t(const hhg::Point3D&)> eight = [](const hhg::Point3D&) { return 8; };
  std::function<real_t(const hhg::Point3D&)> nine = [](const hhg::Point3D&) { return 9; };

  hhg::Face *face0 = storage.beginFaces()->second;

  walberla::mpi::SendBuffer sb;
  hhg::P1BubbleVertex::packData(level,*edge.v0,0,sb);
  hhg::P1BubbleVertex::packData(level,*edge.v1,0,sb);
  walberla::mpi::RecvBuffer rb(sb);
  unpackVertexData(level,edge,memory_id,rb,*edge.v0);
  
  hhg::P1BubbleFace::interpolate(*face0,0,eight,maxLevel);

  WALBERLA_UNUSED(face0);

}