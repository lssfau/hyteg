#include <core/mpi/Environment.h>
#include "tinyhhg_core/tinyhhg.hpp"
#include "tinyhhg_core/p1functionspace/P1DataHandling.hpp"

using namespace hhg;

int main (int argc, char** argv) {
  walberla::mpi::Environment MPIenv( argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  const uint_t minLevel = 3;
  const uint_t maxLevel = 3;

  size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);

  hhg::MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, uint_c(walberla::MPIManager::instance()->numProcesses()));
  RoundRobin loadBalancer;
  setupStorage.balanceLoad( loadBalancer, 0.0 );

  WALBERLA_MPI_BARRIER();

  std::shared_ptr<hhg::PrimitiveStorage> storage(new hhg::PrimitiveStorage(setupStorage));


  auto faceP1FunctionMemoryDataHandling = std::make_shared< FaceP1FunctionMemoryDataHandling >(minLevel, maxLevel);
  auto edgeP1FunctionMemoryDataHandling = std::make_shared< EdgeP1FunctionMemoryDataHandling >(minLevel, maxLevel);
  auto vertexP1FunctionMemoryDataHandling = std::make_shared< VertexP1FunctionMemoryDataHandling >(minLevel, maxLevel);

  PrimitiveDataID<FaceP1FunctionMemory  , Face  > faceDataID;
  PrimitiveDataID<EdgeP1FunctionMemory  , Edge  > edgeDataID;
  PrimitiveDataID<VertexP1FunctionMemory, Vertex> vertexDataID;

  storage->addFaceData(faceDataID, faceP1FunctionMemoryDataHandling, "data");
  storage->addEdgeData(edgeDataID, edgeP1FunctionMemoryDataHandling, "data");
  storage->addVertexData(vertexDataID, vertexP1FunctionMemoryDataHandling, "data");

  hhg::PrimitiveStorage::PrimitiveMap primitveMap;
  storage->getPrimitives( primitveMap );

  std::shared_ptr< communication::BufferedCommunicator > comm = std::make_shared<communication::BufferedCommunicator>( storage ) ;

  WALBERLA_ABORT("TODO");

  //FaceP1BubbleFunctionMemory *face0data = storage->beginFaces()->second->getData(faceDataID);
//  real_t *face0data = storage.get()->beginFaces()->second->getData(faceDataID)->data[maxLevel].get();
//  std::cout << "First Face ID is: " << storage.get()->beginFaces()->second->getID().getID() << std::endl;
//  real_t *face1data = (++storage.get()->beginFaces())->second->getData(faceDataID)->data[maxLevel].get();
//
//  for(auto it = storage->beginEdges(); it != storage->endEdges(); ++it){
//    for(uint_t i = 0; i < v_perEdge; ++i){
//      it->second->getData(edgeDataID)->data[maxLevel].get()[i] = 2. + (walberla::real_c(it->first) * 0.1);
//    }
//
//    if(it->second->getNumHigherDimNeighbors() == 2){
//      std::cout << "Middle Edge ID is: " << it->second->getID().getID() << std::endl;
//    }
//  }
//  for(uint_t i = 0; i < hhg::levelinfo::num_microvertices_per_face(maxLevel);++i){
//    face0data[i] = 1.1;
//    face1data[i] = 1.2;
//  }
////  for(uint_t i = 0; i < hhg::levelinfo::num_microvertices_per_face(maxLevel);++i){
////    face0data[i] = walberla::MPIManager::instance()->rank();
////  }
//  std::shared_ptr<P1BubblePackInfo> packInfo(new P1BubblePackInfo(maxLevel,vertexDataID,edgeDataID,faceDataID,storage));
//  communication::BufferedCommunicator communicator( storage );
//  communicator.addPackInfo(packInfo);
//  communicator.startCommunication<Edge,Face>();
//  communicator.endCommunication<Edge,Face>();
//
//
//  std::function<real_t(const hhg::Point3D&)> eight = [](const hhg::Point3D&) { return 8; };
//  std::function<real_t(const hhg::Point3D&)> nine = [](const hhg::Point3D&) { return 9; };
//
////  for(auto it = storage->beginEdges(); it != storage->endEdges(); ++it){
////    walberla::mpi::SendBuffer sb;
////    for(auto jt = it->second->beginHigherDimNeighbors(); jt != it->second->endHigherDimNeighbors(); ++jt){
////      packInfo->packEdgeForFace(it->second,jt->first,sb);
////      walberla::mpi::RecvBuffer rb(sb);
////      packInfo->unpackFaceFromEdge(storage->getFace(jt->first),it->first,rb);
////    }
////  }
//
//  //hhg::P1BubbleFace::interpolate(*face0,0,eight,maxLevel);
//
//
//
//  std::cout << "++++ Face 0 Vertex Dofs: ++++" << std::endl;
//  for(size_t i = 0; i < v_perEdge; ++i){
//    for(size_t j = 0; j < v_perEdge - i; ++j) {
//      //std::cout << face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
//      fmt::print("{0:.2f}  ",face0data[hhg::P1BubbleFace::CoordsVertex::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsVertex::VERTEX_C)]);
//      //std::cout << face0data[hhg::P1BubbleFace::CoordsVertex::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsVertex::VERTEX_C)];
//    }
//    std::cout << std::endl;
//  }
//  std::cout << "++++ Face 1 Vertex Dofs: ++++" << std::endl;
//  for(size_t i = 0; i < v_perEdge; ++i){
//    for(size_t j = 0; j < v_perEdge - i; ++j) {
//      //std::cout << face0mem[CoordsVertex::index<maxLevel>(i, j, CoordsVertex::VERTEX_C)] << " ";
//      fmt::print("{0:.2f}  ",face1data[hhg::P1BubbleFace::CoordsVertex::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsVertex::VERTEX_C)]);
//      //std::cout << face0data[hhg::P1BubbleFace::CoordsVertex::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsVertex::VERTEX_C)];
//    }
//    std::cout << std::endl;
//  }
}
