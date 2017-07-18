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

  size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);

  hhg::MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, 1  );

  std::shared_ptr<hhg::PrimitiveStorage> storage(new hhg::PrimitiveStorage(uint_c(walberla::MPIManager::instance()->rank()) , setupStorage));

  FaceP1BubbleFunctionMemoryDataHandling faceP1BubbleFunctionMemoryDataHandling;
  EdgeP1BubbleFunctionMemoryDataHandling edgeP1BubbleFunctionMemoryDataHandling;
  VertexP1BubbleFunctionMemoryDataHandling vertexP1BubbleFunctionMemoryDataHandling;
  hhg::PrimitiveDataID<FaceP1BubbleFunctionMemory  , Face  > faceDataID   = storage->addFaceData(faceP1BubbleFunctionMemoryDataHandling,"data");
  hhg::PrimitiveDataID<EdgeP1BubbleFunctionMemory  , Edge  > edgeDataID   = storage->addEdgeData(edgeP1BubbleFunctionMemoryDataHandling,"data");
  hhg::PrimitiveDataID<VertexP1BubbleFunctionMemory, Vertex> vertexDataID = storage->addVertexData(vertexP1BubbleFunctionMemoryDataHandling,"data");

  hhg::PrimitiveStorage::PrimitiveMap primitveMap;
  storage->getPrimitives( primitveMap);

  //FaceP1BubbleFunctionMemory *face0data = storage->beginFaces()->second->getData(faceDataID);
  real_t *face0data = storage.get()->beginFaces()->second->getData(faceDataID)->data[maxLevel].get();
  std::cout << "First Face ID is: " << storage.get()->beginFaces()->second->getID().getID() << std::endl;
  real_t *face1data = (++storage.get()->beginFaces())->second->getData(faceDataID)->data[maxLevel].get();

  for(auto it = storage->beginEdges(); it != storage->endEdges(); ++it){
    for(uint_t i = 0; i < v_perEdge; ++i){
      it->second->getData(edgeDataID)->data[maxLevel].get()[i] = 2 + (it->first * 0.1);
    }

    if(it->second->getNumHigherDimNeighbors() == 2){
      EdgeP1BubbleFunctionMemory *middleEdgeData = it->second->getData(edgeDataID);
      std::cout << "Middle Edge ID is: " << it->second->getID().getID() << std::endl;
    }
  }
  real_t counter = 1;
  for(uint_t i = 0; i < hhg::levelinfo::num_microvertices_per_face(maxLevel);++i){
    face0data[i] = 1.1;
  }
  P1BubblePackInfo *packInfo = new P1BubblePackInfo(maxLevel,vertexDataID,edgeDataID,faceDataID,storage);
//  for(auto it = storage->beginEdges(); it != storage->endEdges(); ++it){
//
//    //packInfo->packEdgeForFace(*it,)
//  }


  std::function<real_t(const hhg::Point3D&)> eight = [](const hhg::Point3D&) { return 8; };
  std::function<real_t(const hhg::Point3D&)> nine = [](const hhg::Point3D&) { return 9; };
//
//  hhg::Face *face0 = storage->beginFaces()->second;
//
//  walberla::mpi::SendBuffer sb;
//  packInfo->packEdgeForFace(storage->getEdge(7),10, sb);
//  hhg::P1BubbleVertex::packData(level,*edge.v0,0,sb);
// hhg::P1BubbleVertex::packData(level,*edge.v1,0,sb);
//  walberla::mpi::RecvBuffer rb(sb);
//  packInfo->unpackFaceFromEdge(storage->getFace(10),7,rb);

  for(auto it = storage->beginEdges(); it != storage->endEdges(); ++it){
    walberla::mpi::SendBuffer sb;
    for ( const auto & higherDimNeighborID : it->second->higherDimNeighbors() ) {
      packInfo->packEdgeForFace(it->second,higherDimNeighborID,sb);
      walberla::mpi::RecvBuffer rb(sb);
      packInfo->unpackFaceFromEdge(storage->getFace(higherDimNeighborID),it->first,rb);
    }
  }

  //hhg::P1BubbleFace::interpolate(*face0,0,eight,maxLevel);



  std::cout << "++++ Face 0 Vertex Dofs: ++++" << std::endl;
  for(size_t i = 0; i < v_perEdge; ++i){
    for(size_t j = 0; j < v_perEdge - i; ++j) {
      //std::cout << face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
      fmt::print("{0:.2f}  ",face0data[hhg::P1BubbleFace::CoordsVertex::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsVertex::VERTEX_C)]);
      //std::cout << face0data[hhg::P1BubbleFace::CoordsVertex::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsVertex::VERTEX_C)];
    }
    std::cout << std::endl;
  }
  std::cout << "++++ Face 1 Vertex Dofs: ++++" << std::endl;
  for(size_t i = 0; i < v_perEdge; ++i){
    for(size_t j = 0; j < v_perEdge - i; ++j) {
      //std::cout << face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
      fmt::print("{0:.2f}  ",face1data[hhg::P1BubbleFace::CoordsVertex::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsVertex::VERTEX_C)]);
      //std::cout << face0data[hhg::P1BubbleFace::CoordsVertex::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsVertex::VERTEX_C)];
    }
    std::cout << std::endl;
  }
}
