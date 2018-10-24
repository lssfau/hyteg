#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilDirections.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/indexing/MacroFaceIndexing.hpp"


#include "core/mpi/all.h"
#include "core/debug/all.h"

using namespace hhg;

using walberla::real_t;

void checkComm3d( const uint_t level )
{
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   hhg::EdgeDoFFunction< int > x( "x", storage, level, level );
   /// for y we set the local comm to mpi; default would be direct
   hhg::EdgeDoFFunction< int > y( "x", storage, level, level );
   y.setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );

   ////////// check cell to face comm //////////
   for( auto& cellIt : storage->getCells() )
   {
      Cell& cell      = *cellIt.second;
      int*  cellData  = cell.getData( x.getCellDataID() )->getPointer( level );
      int*  cellDataY = cell.getData( y.getCellDataID() )->getPointer( level );
      for( uint_t i = 0; i < cell.getData( x.getCellDataID() )->getSize( level ); ++i )
      {
         cellData[i]  = 13;
         cellDataY[i] = 26;
      }
   }

   x.communicate< Cell, Face >( level );
   y.communicate< Cell, Face >( level );

   for( auto& faceIt : storage->getFaces() )
   {
      Face& face      = *faceIt.second;
      int*  faceData  = face.getData( x.getFaceDataID() )->getPointer( level );
      int*  faceDataY = face.getData( y.getFaceDataID() )->getPointer( level );
      /// all non inner DoFs should be set to 13/26 so we start after the inner DoFs
      for( uint_t i = hhg::levelinfo::num_microedges_per_face( level ); i < face.getData( x.getFaceDataID() )->getSize( level );
           ++i )
      {
         WALBERLA_CHECK_EQUAL( faceData[i], 13, i );
         WALBERLA_CHECK_EQUAL( faceDataY[i], 26, i );
      }
   }
   /////////////////////////////////////////////

   ////////// check edge to face comm //////////
   for( auto& edgeIt : storage->getEdges() )
   {
      Edge& edge      = *edgeIt.second;
      int*  edgeData  = edge.getData( x.getEdgeDataID() )->getPointer( level );
      int*  edgeDataY = edge.getData( y.getEdgeDataID() )->getPointer( level );
      for( uint_t i = 0; i < levelinfo::num_microedges_per_edge( level ); ++i )
      {
         edgeData[i]  = 14;
         edgeDataY[i] = 26;
      }
   }

   x.communicate< Edge, Face >( level );
   y.communicate< Edge, Face >( level );

   using hhg::edgedof::macroface::BorderIterator;

   for( auto& faceIt : storage->getFaces() )
   {
      Face& face      = *faceIt.second;
      int*  faceData  = face.getData( x.getFaceDataID() )->getPointer( level );
      int*  faceDataY = face.getData( y.getFaceDataID() )->getPointer( level );
      for( const auto& it : BorderIterator( level , indexing::FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT, 0 ) )
      {
         WALBERLA_CHECK_EQUAL( faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_HO_C )], 14);
         WALBERLA_CHECK_EQUAL( faceDataY[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_HO_C )], 26);
      }

      for( const auto& it : BorderIterator( level , indexing::FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 0 ) )
      {
         WALBERLA_CHECK_EQUAL( faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_DI_N )], 14);
         WALBERLA_CHECK_EQUAL( faceDataY[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_DI_N )], 26);
      }

      for( const auto& it : BorderIterator( level , indexing::FaceBorderDirection::LEFT_BOTTOM_TO_TOP, 0 ) )
      {
         WALBERLA_CHECK_EQUAL( faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_VE_NW )], 14);
         WALBERLA_CHECK_EQUAL( faceDataY[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_VE_NW )], 26);
      }
   }

   /////////////////////////////////////////////

//   ////////// check face to edge comm //////////
//   for( auto& faceIt : storage->getFaces() )
//   {
//      Face& face      = *faceIt.second;
//      int*  faceData  = face.getData( x.getFaceDataID() )->getPointer( level );
//      int*  faceDataY = face.getData( y.getFaceDataID() )->getPointer( level );
//      for( uint_t i = 0; i < face.getData( x.getFaceDataID() )->getSize( level ); ++i )
//      {
//         faceData[i]  = 17;
//         faceDataY[i] = 33;
//      }
//   }
//
//   x.communicate< Face, Edge >( level );
//   y.communicate< Face, Edge >( level );
//
//   for( auto& edgeIt : storage->getEdges() )
//   {
//      Edge& edge      = *edgeIt.second;
//      int*  edgeData  = edge.getData( x.getEdgeDataID() )->getPointer( level );
//      int*  edgeDataY = edge.getData( y.getEdgeDataID() )->getPointer( level );
//      for( uint_t i = 0; i < levelinfo::num_microedges_per_edge( level ); ++i )
//      {
//         WALBERLA_CHECK_EQUAL( edgeData[i], 17, i );
//         WALBERLA_CHECK_EQUAL( edgeDataY[i], 33, i );
//      }
//   }
//   /////////////////////////////////////////////
}

template<uint_t level>
void checkComm(const std::string &meshfile, bool bufferComm = false){

  MeshInfo meshInfo = MeshInfo::fromGmshFile(meshfile);
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  hhg::EdgeDoFFunction< int > x("x", storage, level, level);
  if(bufferComm) {
    x.setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
  }

  x.enumerate(level);

  uint_t numberOfChecks = 0;
  uint_t totalExpectedChecks = 0;

  for(auto &edgeIt : storage->getEdges()){
    if(edgeIt.second.get()->getNumHigherDimNeighbors() == 1){
      totalExpectedChecks += 3 * levelinfo::num_microedges_per_edge( level ) + levelinfo::num_microedges_per_edge( level ) - 1;
    } else if(edgeIt.second.get()->getNumHigherDimNeighbors() == 2){
      totalExpectedChecks += 6 * levelinfo::num_microedges_per_edge( level ) + 2 * (levelinfo::num_microedges_per_edge( level ) - 1);
    } else {
      WALBERLA_CHECK(false);
    }
  }
  for(auto &vertexIt : storage->getVertices()){
    totalExpectedChecks += vertexIt.second->getNumNeighborFaces();
    totalExpectedChecks += vertexIt.second->getNumNeighborEdges();
  }

  using hhg::edgedof::macroface::BorderIterator;
  for (auto &faceIt : storage->getFaces()) {
    Face &face = *faceIt.second;
    int *faceData = face.getData(x.getFaceDataID())->getPointer(level);
    std::vector<PrimitiveID> nbrEdges;
    face.getNeighborEdges(nbrEdges);
    uint_t localEdgeIdOnFace = 0;

/////////// FIRST EDGE ////////////

    Edge *firstEdge = storage->getEdge(nbrEdges[localEdgeIdOnFace].getID());
    int *edgeData = firstEdge->getData(x.getEdgeDataID())->getPointer(level);
    uint_t idxCounter = 0;
    /// horizontal Dof on edge 0
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),0)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromHorizontalEdge( level, idxCounter, stencilDirection::EDGE_HO_C )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_HO_C )]
      , "it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter)
      idxCounter++;
      numberOfChecks++;
    }
    /// horizontal Dof on Face for edge 0; offset 1 to border
    idxCounter = 1;
    stencilDirection edgeDir = firstEdge->face_index(face.getID()) == 0 ? stencilDirection::EDGE_HO_SE : stencilDirection::EDGE_HO_NW;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),1)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromVertex( level, idxCounter, edgeDir )],
              faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_HO_C )]
            ,"it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter);
      idxCounter++;
      numberOfChecks++;

    }
    /// diagonal Dof on Face for edge 0; offset 0 to border
    idxCounter = 1;
    edgeDir = firstEdge->face_index(face.getID()) == 0 ? stencilDirection::EDGE_VE_S : stencilDirection::EDGE_DI_NW;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),0)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromVertex( level, idxCounter, edgeDir )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_DI_N )]
      ,"it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter);
      idxCounter++;
      numberOfChecks++;
    }
    /// vertical Dof on Face for edge 0; offset 0 to border
    idxCounter = 1;
    edgeDir = firstEdge->face_index(face.getID()) == 0 ? stencilDirection::EDGE_DI_SW : stencilDirection::EDGE_VE_NW;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),0)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromVertex( level, idxCounter, edgeDir )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_VE_NW )]
      ,"it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter);
      idxCounter++;
      numberOfChecks++;
    }
/////////// SECOND EDGE ////////////
    localEdgeIdOnFace = 1;
    Edge *secondEdge = storage->getEdge(nbrEdges[localEdgeIdOnFace].getID());
    edgeData = secondEdge->getData(x.getEdgeDataID())->getPointer(level);
    /// horizontal Dof on edge 1
    idxCounter = 0;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),0)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromHorizontalEdge( level, idxCounter, stencilDirection::EDGE_HO_C )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_DI_N )]
      , "it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter)
      idxCounter++;
      numberOfChecks++;
    }
    /// diagonal Dof on Face = horizontal Dof on edge; offset 1 to border
    idxCounter = 1;
    edgeDir = secondEdge->face_index(face.getID()) == 0 ? stencilDirection::EDGE_HO_SE : stencilDirection::EDGE_HO_NW;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),1)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromVertex( level, idxCounter, edgeDir )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_DI_N )]
      ,"it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter);
      idxCounter++;
      numberOfChecks++;
    }
    /// vertical Dof on Face = diagonal Dof on edge; offset 1 to border
    idxCounter = 1;
    edgeDir = secondEdge->face_index(face.getID()) == 0 ? stencilDirection::EDGE_VE_S : stencilDirection::EDGE_DI_NW;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),0)) {
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromVertex( level, idxCounter, edgeDir )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_VE_NW )],
        "it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter);
      idxCounter++;
      numberOfChecks++;
    }
      /// horizontal Dof on Face = vertical Dof on edge; offset 1 to border
      idxCounter = 1;
      edgeDir = secondEdge->face_index(face.getID()) == 0 ? stencilDirection::EDGE_DI_SW : stencilDirection::EDGE_VE_NW;
      for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),0)){
        WALBERLA_CHECK_EQUAL(
          edgeData[edgedof::macroedge::indexFromVertex( level, idxCounter, edgeDir )],
          faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_HO_C )]
        ,"it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter);
        idxCounter++;
        numberOfChecks++;
      }
/////////// THIRD EDGE ////////////
    localEdgeIdOnFace = 2;
    Edge *thirdEdge = storage->getEdge(nbrEdges[localEdgeIdOnFace].getID());
    edgeData = thirdEdge->getData(x.getEdgeDataID())->getPointer(level);
    /// horizontal Dof on edge 2
    idxCounter = 0;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),0)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromHorizontalEdge( level, idxCounter, stencilDirection::EDGE_HO_C )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_VE_NW )]
      , "it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter)
      idxCounter++;
      numberOfChecks++;
    }
    /// vertical Dof on face for edge 2 = horizontal on edge; offset 1 to border
    idxCounter = 1;
    edgeDir = thirdEdge->face_index(face.getID()) == 0 ? stencilDirection::EDGE_HO_SE : stencilDirection::EDGE_HO_NW;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),1)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromVertex( level, idxCounter, edgeDir )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_VE_NW )]
      ,"it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter);
      idxCounter++;
      numberOfChecks++;
    }
    /// horizontal Dof on face for edge 2 = diagonal on edge; offset 1 to border
    idxCounter = 1;
    edgeDir = thirdEdge->face_index(face.getID()) == 0 ? stencilDirection::EDGE_DI_SW : stencilDirection::EDGE_VE_NW;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),0)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromVertex( level, idxCounter, edgeDir )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_HO_C )]
      ,"it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter);
      idxCounter++;
      numberOfChecks++;
    }
    /// diagonal Dof on face for edge 2 = vertical on edge; offset 1 to border
    idxCounter = 1;
    edgeDir = thirdEdge->face_index(face.getID()) == 0 ? stencilDirection::EDGE_VE_S : stencilDirection::EDGE_DI_NW;
    for(const auto& it : BorderIterator(level,indexing::getFaceBorderDirection(localEdgeIdOnFace,face.edge_orientation[localEdgeIdOnFace]),0)){
      WALBERLA_CHECK_EQUAL(
        edgeData[edgedof::macroedge::indexFromVertex( level, idxCounter, edgeDir )],
        faceData[edgedof::macroface::indexFromHorizontalEdge( level, it.col(), it.row(), stencilDirection::EDGE_DI_N )]
      ,"it.col(): " << it.col() << " it.row(): " << it.row() << " idxCounter: " << idxCounter);
      idxCounter++;
      numberOfChecks++;
    }
  }

  for (auto &vertexIt : storage->getVertices()) {
    Vertex &vertex = *vertexIt.second;
    int *vertexData = vertex.getData(x.getVertexDataID())->getPointer(level);

    for (const PrimitiveID& edgeId : vertex.neighborEdges()) {
      Edge *edge = storage->getEdge(edgeId);
      int *edgeData = edge->getData(x.getEdgeDataID())->getPointer(level);
      if (edge->getVertexID0() == vertex.getID()) {
        WALBERLA_CHECK_EQUAL(
          edgeData[edgedof::macroedge::indexFromVertex( level, 1, stencilDirection::EDGE_HO_W )],
          vertexData[vertex.edge_index(edgeId)],
          "vertex: " << vertex.getID().getID() << " edgeIndex: " << vertex.edge_index(edgeId))
        numberOfChecks++;
      } else if (edge->getVertexID1() == vertex.getID()) {
        WALBERLA_CHECK_EQUAL(
          edgeData[edgedof::macroedge::indexFromVertex( level, levelinfo::num_microvertices_per_edge( level ) - 1, stencilDirection::EDGE_HO_W )],
          vertexData[vertex.edge_index(edgeId)],
          " edgeIndex: " << vertex.edge_index(edgeId))
        numberOfChecks++;
      } else {
        WALBERLA_ABORT("edge is not on vertex")
      }
    }
    for (const PrimitiveID& faceId : vertex.neighborFaces()) {
      Face *face = storage->getFace(faceId);
      int *faceData = face->getData(x.getFaceDataID())->getPointer(level);
      if (face->getVertexID0() == vertex.getID()) {
        WALBERLA_CHECK_EQUAL(
          faceData[edgedof::macroface::indexFromDiagonalEdge( level, 0, 0, stencilDirection::EDGE_DI_C )],
          vertexData[vertex.getNumNeighborEdges() + vertex.face_index(faceId)],
          " faceIndex: " << vertex.face_index(faceId))
        numberOfChecks++;
      } else if (face->getVertexID1() == vertex.getID()){
        uint_t nbrEdgeDoFs = levelinfo::num_microedges_per_edge( level );
        WALBERLA_CHECK_EQUAL(
          faceData[edgedof::macroface::indexFromVerticalEdge( level, nbrEdgeDoFs - 1, 0, stencilDirection::EDGE_VE_C )],
          vertexData[vertex.getNumNeighborEdges() + vertex.face_index(faceId)],
          " index: " << vertex.getNumNeighborEdges() + vertex.face_index(faceId))
        numberOfChecks++;
      } else if (face->getVertexID2() == vertex.getID()){
        uint_t nbrEdgeDoFs = levelinfo::num_microedges_per_edge( level );
        WALBERLA_CHECK_EQUAL(
          faceData[edgedof::macroface::indexFromHorizontalEdge( level, 0, nbrEdgeDoFs - 1, stencilDirection::EDGE_HO_C )],
          vertexData[vertex.getNumNeighborEdges() + vertex.face_index(faceId)],
          " faceIndex: " << vertex.getNumNeighborEdges() + vertex.face_index(faceId))
        numberOfChecks++;
      } else {
        WALBERLA_ABORT("face it not on vertex");
      }
    }

  }

  WALBERLA_CHECK_EQUAL(totalExpectedChecks,numberOfChecks, "expected: " << totalExpectedChecks << " number: " << numberOfChecks);

}

int main (int argc, char ** argv ) {

  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  walberla::debug::enterTestMode();

//  checkComm("../../data/meshes/tri_1el.msh", true);
//
//  checkComm("../../data/meshes/tri_1el.msh", false);
  checkComm<3>("../../data/meshes/tri_1el.msh", true);

  checkComm<3>("../../data/meshes/tri_1el.msh", false);

  checkComm<4>("../../data/meshes/tri_1el.msh", true);

  checkComm<4>("../../data/meshes/tri_1el.msh", false);


  checkComm<3>("../../data/meshes/quad_4el.msh", true);

  checkComm<4>("../../data/meshes/quad_4el.msh", true);

  checkComm<5>("../../data/meshes/quad_4el.msh", true);

  checkComm<4>("../../data/meshes/quad_4el.msh", false);

  checkComm<5>("../../data/meshes/quad_4el.msh", false);

  checkComm<3>("../../data/meshes/bfs_12el.msh", true);

  checkComm<3>("../../data/meshes/bfs_12el.msh", false);

  checkComm3d( 2u );

   return 0;


}
