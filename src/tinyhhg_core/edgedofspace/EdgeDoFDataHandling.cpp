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


}/// namespace hhg