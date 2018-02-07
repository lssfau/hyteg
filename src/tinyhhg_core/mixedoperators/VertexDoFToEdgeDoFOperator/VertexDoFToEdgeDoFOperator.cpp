#include "VertexDoFToEdgeDoFOperator.hpp"

#include "tinyhhg_core/p2functionspace/P2Elements.hpp"

namespace hhg {

template<class UFCOperator>
VertexDoFToEdgeDoFOperator<UFCOperator>::VertexDoFToEdgeDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel, size_t maxLevel)
  : Operator(storage, minLevel, maxLevel) {
  /// since the Vertex does not own any EdgeDoFs only edge and face are needed

  auto edgeDataHandling =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Edge >>(minLevel_,
                                                                        maxLevel_,
                                                                        VertexDoFToEdgeDoF::macroEdgeVertexDoFToEdgeDoFStencilSize);

  auto faceDataHandling =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Face >>(minLevel_,
                                                                        maxLevel_,
                                                                        VertexDoFToEdgeDoF::macroFaceVertexDoFToEdgeDoFStencilSize);

  storage->addEdgeData(edgeStencilID_, edgeDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addFaceData(faceStencilID_, faceDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil");

  // Only assemble stencils if UFCOperator is specified
  if (!std::is_same<UFCOperator, fenics::NoAssemble>::value) {
    assembleStencils();
  }
}

template<class UFCOperator>
void VertexDoFToEdgeDoFOperator<UFCOperator>::assembleStencils() {
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

      // Assemble vertexToEdge stencil
      real_t * vStencil = storage_->getFace(face.getID())->getData(faceStencilID_)->getPointer(level);
      P2Face::VertexToEdge::assembleStencil(local_stiffness_gray, local_stiffness_blue, vStencil);
//      WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Face = {}", PointND<real_t, 12>(&vStencil[0])));
    }

    // Assemble edge stencils
    for (auto& it : storage_->getEdges()) {
      Edge &edge = *it.second;

      // Assemble vertexToEdge stencil
      Face* face = storage_->getFace(edge.neighborFaces()[0]);
      real_t* vStencil = storage_->getEdge(edge.getID())->getData(edgeStencilID_)->getPointer(level);
      compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
      compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
      P2Edge::VertexToEdge::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true);

      if (edge.getNumNeighborFaces() == 2) {
        face = storage_->getFace(edge.neighborFaces()[1]);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
        P2Edge::VertexToEdge::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false);
      }

//      WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Edge = {}", PointND<real_t, 4>(&vStencil[0])));
    }

  }
}

template<class UFCOperator>
void VertexDoFToEdgeDoFOperator<UFCOperator>::compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type) {
  real_t coords[6];
  fenics::compute_micro_coords(face, level, coords, element_type);
  UFCOperator gen;
  gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
}

template<class UFCOperator>
void VertexDoFToEdgeDoFOperator<UFCOperator>::apply_impl(P1Function<real_t> &src, EdgeDoFFunction<real_t> &dst, size_t level, DoFType flag,
                                            UpdateType updateType) {

  ///the order of communication is crucial here.
  ///first the vertex dofs on the macro vertex need to be communicated to the edge since they are needed on the edge and the face
  src.getCommunicator(level)->communicate<Vertex, Edge>();
  ///secondly the vertex dofs on the macro edge are communicated to the face passing on the vertex dof from the macro vertex
  src.getCommunicator(level)->communicate<Edge, Face>();
  ///lastly the vertex dofs on the macro face are communicated to the edge which also contain vertex dofs which are located on neighboring edges
  src.getCommunicator(level)->startCommunication<Face, Edge>();

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      VertexDoFToEdgeDoF::applyFace(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
    }
  }

  src.getCommunicator(level)->endCommunication<Face, Edge>();

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      VertexDoFToEdgeDoF::applyEdge(level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
    }
  }

}

namespace VertexDoFToEdgeDoF {

///
/// \param level stencil size is independent of level
/// \param numDependencies number of adjacent faces of the edge
/// \return number of the stencil entries
uint_t macroEdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies)
{
  WALBERLA_UNUSED( level );
  return 2 + numDependencies;
}

uint_t macroFaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies)
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( numDependencies );
  return 4 + 4 + 4;
}
}

template class VertexDoFToEdgeDoFOperator<hhg::fenics::NoAssemble>;
template class VertexDoFToEdgeDoFOperator<p2_divt_cell_integral_0_otherwise>;
template class VertexDoFToEdgeDoFOperator<p2_divt_cell_integral_1_otherwise>;

}/// namespace hhg