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


/// on edges only one stencil is required since only the horizontal edge DoFs belong to the edge
uint_t macroEdgeEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 1 + 2 * numDependencies;
}
/// on face three stencils are needed for horizontal, vertical and diagonal DoFs
uint_t macroFaceEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  WALBERLA_UNUSED(numDependencies);
  return 5 + 5 + 5;
}

////////////////////
// Implementation //
////////////////////


std::shared_ptr<StencilMemory< real_t >> MacroEdgeEdgeDoFToEdgeDoFDataHandling::initialize(const Edge *const edge) const {
  return std::make_shared<StencilMemory<real_t> >(macroEdgeEdgeDoFToEdgeDoFStencilSize,
                                                  edge->getNumNeighborFaces(),
                                                  minLevel_,
                                                  maxLevel_);
}


std::shared_ptr<StencilMemory< real_t >> MacroFaceEdgeDoFToEdgeDoFDataHandling::initialize(const Face *const face) const {
  WALBERLA_UNUSED(face);
  return std::make_shared<StencilMemory<real_t> >(macroFaceEdgeDoFToEdgeDoFStencilSize,
                                                  0,
                                                  minLevel_,
                                                  maxLevel_);
}


}/// namespace hhg