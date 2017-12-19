#include "EdgeDoFDataHandling.hpp"

namespace hhg {

//////////////////////
/// Size Functions ///
//////////////////////

uint_t edgeDoFMacroVertexFunctionMemorySize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 2 * numDependencies;
}

uint_t edgeDoFMacroEdgeFunctionMemorySize(const uint_t &level, const uint_t &numDependencies) {
  return levelinfo::num_microedges_per_edge(level) + numDependencies * (3 * (levelinfo::num_microedges_per_edge(level)) - 1);
}

uint_t edgeDoFMacroFaceFunctionMemorySize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(numDependencies);
  return 3 * (((levelinfo::num_microedges_per_edge(level) + 1) * levelinfo::num_microedges_per_edge(level)) / 2);
}

uint_t macroVertexEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 2 * numDependencies;
}

uint_t macroEdgeEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 2 + 5 * numDependencies;
}

uint_t macroFaceEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  WALBERLA_UNUSED(numDependencies);
  return 12;
}

////////////////////
// Implementation //
////////////////////

std::shared_ptr<StencilMemory<double>> MacroVertexEdgeDoFToEdgeDoFDataHandling::initialize(const Vertex *const vertex) const {
  return std::make_shared<StencilMemory< real_t > >( macroVertexEdgeDoFToEdgeDoFStencilSize,
                                                     vertex->getNumNeighborEdges(),
                                                     minLevel_,
                                                     maxLevel_);
}

std::shared_ptr<StencilMemory<double>> MacroEdgeEdgeDoFToEdgeDoFDataHandling::initialize(const Edge *const edge) const {
  return std::make_shared<StencilMemory<real_t> >(macroEdgeEdgeDoFToEdgeDoFStencilSize,
                                                  edge->getNumNeighborFaces(),
                                                  minLevel_,
                                                  maxLevel_);
}


std::shared_ptr<StencilMemory<double>> MacroFaceEdgeDoFToEdgeDoFDataHandling::initialize(const Face *const face) const {
  return std::make_shared<StencilMemory<real_t> >(macroFaceEdgeDoFToEdgeDoFStencilSize,
                                                  0,
                                                  minLevel_,
                                                  maxLevel_);
}


}/// namespace hhg