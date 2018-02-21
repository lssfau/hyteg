#include "EdgeDoFToVertexDoFOperator.hpp"
#include "EdgeDoFToVertexDoFApply.hpp"

#include "tinyhhg_core/p2functionspace/P2Elements.hpp"

namespace hhg{

template<class UFCOperator>
EdgeDoFToVertexDoFOperator<UFCOperator>::EdgeDoFToVertexDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
                                                       const size_t & minLevel,
                                                       const size_t & maxLevel)
  :Operator(storage,minLevel,maxLevel)
{

  using namespace EdgeDoFToVertexDoF;

  auto vertexDataHandling =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Vertex >>(minLevel_, maxLevel_, macroVertexEdgeDoFToVertexDoFStencilSize);

  auto edgeDataHandling   =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Edge   >>(minLevel_, maxLevel_, macroEdgeEdgeDoFToVertexDoFStencilSize);

  auto faceDataHandling   =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Face   >>(minLevel_, maxLevel_, macroFaceEdgeDoFToVertexDoFStencilSize);

  storage->addVertexData(vertexStencilID_, vertexDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addEdgeData(edgeStencilID_, edgeDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
  storage->addFaceData(faceStencilID_, faceDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");

  // Only assemble stencils if UFCOperator is specified
  if (!std::is_same<UFCOperator, fenics::NoAssemble>::value) {
    assembleStencils();
  }
}

template<class UFCOperator>
void EdgeDoFToVertexDoFOperator<UFCOperator>::assembleStencils() {
  using namespace P2Elements;

  // Initialize memory for local 6x6 matrices
  Matrix6r local_stiffness_gray;
  Matrix6r local_stiffness_blue;

  // Assemble stencils on all levels
  for (uint_t level = minLevel_; level <= maxLevel_; ++level)
  {

    // Assemble face stencils
    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      // Compute both local stiffness matrices
      compute_local_stiffness(face, level, local_stiffness_gray, fenics::GRAY);
      compute_local_stiffness(face, level, local_stiffness_blue, fenics::BLUE);

      // Assemble edgeToVertex stencil
      real_t* vStencil = storage_->getFace(face.getID())->getData(faceStencilID_)->getPointer(level);
      P2Face::EdgeToVertex::assembleStencil(local_stiffness_gray, local_stiffness_blue, vStencil);
//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Face = {}", PointND<real_t, 12>(&vStencil[0])));
    }

    // Assemble edge stencils
    for (auto& it : storage_->getEdges()) {
      Edge &edge = *it.second;

      // Assemble edgeToVertex
      Face* face = storage_->getFace(edge.neighborFaces()[0]);
      real_t* vStencil = storage_->getEdge(edge.getID())->getData(edgeStencilID_)->getPointer(level);
      compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
      compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
      P2Edge::EdgeToVertex::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true);

      if (edge.getNumNeighborFaces() == 2) {
        face = storage_->getFace(edge.neighborFaces()[1]);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
        P2Edge::EdgeToVertex::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false);
      }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Edge = {}", PointND<real_t, 7>(&vStencil[0])));
    }

    for (auto& it : storage_->getVertices()) {
      Vertex &vertex = *it.second;

      // Assemble EdgeToVertex
      real_t* vStencil = storage_->getVertex(vertex.getID())->getData(vertexStencilID_)->getPointer(level);
      for (auto& faceId : vertex.neighborFaces())
      {
        Face* face = storage_->getFace(faceId);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        P2Vertex::EdgeToVertex::assembleStencil(vertex, *face, local_stiffness_gray, vStencil, storage_);
      }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Vertex = {}", PointND<real_t, 5>(&vStencil[0])));
    }

  }
}

template<class UFCOperator>
void EdgeDoFToVertexDoFOperator<UFCOperator>::compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type) {
  real_t coords[6];
  fenics::compute_micro_coords(face, level, coords, element_type);
  UFCOperator gen;
  gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
}

template<class UFCOperator>
void EdgeDoFToVertexDoFOperator<UFCOperator>::apply_impl(EdgeDoFFunction<real_t> &src,
                                            P1Function<real_t> &dst,
                                            uint_t level,
                                            DoFType flag,
                                            UpdateType updateType)
{
  using namespace EdgeDoFToVertexDoF;
  ///there might be room for optimization in the communication. i.e. splitting communicate into start and end to overlap comm and calc
  src.getCommunicator(level)->communicate<Edge, Face>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      applyFace(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
    }
  }


  src.getCommunicator(level)->communicate<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      applyEdge(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
    }
  }

  src.getCommunicator(level)->communicate<Edge, Vertex>();

  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      applyVertex(level, vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), updateType);
    }
  }

}

template<class UFCOperator>
const PrimitiveDataID<StencilMemory< real_t >, Vertex > &EdgeDoFToVertexDoFOperator<UFCOperator>::getVertexStencilID() const {
  return vertexStencilID_;
}

template<class UFCOperator>
const PrimitiveDataID<StencilMemory< real_t >, Edge > &EdgeDoFToVertexDoFOperator<UFCOperator>::getEdgeStencilID() const {
  return edgeStencilID_;
}

template<class UFCOperator>
const PrimitiveDataID<StencilMemory< real_t >, Face > &EdgeDoFToVertexDoFOperator<UFCOperator>::getFaceStencilID() const {
  return faceStencilID_;
}


namespace EdgeDoFToVertexDoF {
////////// Stencil sizes //////////
uint_t macroVertexEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive ) {
  WALBERLA_UNUSED(level);
  return 2 * primitive.getNumNeighborEdges();
}

uint_t macroEdgeEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive ) {
  WALBERLA_UNUSED(level);
  return 2 + 5 * primitive.getNumNeighborFaces();
}

uint_t macroFaceEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive ) {
  WALBERLA_UNUSED(level);
  WALBERLA_UNUSED(primitive);
  return 12;
}

}/// namespace EdgeDoFToVertexDoF

template class EdgeDoFToVertexDoFOperator<hhg::fenics::NoAssemble>;
template class EdgeDoFToVertexDoFOperator<p2_div_cell_integral_0_otherwise>;
template class EdgeDoFToVertexDoFOperator<p2_div_cell_integral_1_otherwise>;

}/// namespace hhg
