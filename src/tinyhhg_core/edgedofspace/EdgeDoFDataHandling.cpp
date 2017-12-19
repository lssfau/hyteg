#include "EdgeDoFDataHandling.hpp"

namespace hhg {

//////////////////////
/// Size Functions ///
//////////////////////

uint_t EdgeDoFMacroVertexFunctionMemorySize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 2 * numDependencies;
}

uint_t EdgeDoFMacroEdgeFunctionMemorySize(const uint_t &level, const uint_t &numDependencies) {
  return levelinfo::num_microedges_per_edge(level) + numDependencies * (3 * (levelinfo::num_microedges_per_edge(level)) - 1);
}

uint_t EdgeDoFMacroFaceFunctionMemorySize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(numDependencies);
  return 3 * (((levelinfo::num_microedges_per_edge(level) + 1) * levelinfo::num_microedges_per_edge(level)) / 2);
}

////////////////////
// Implementation //
////////////////////

template<typename ValueType>
std::shared_ptr<FunctionMemory<ValueType> >
EdgeDoFMacroVertexFunctionMemoryDataHandling<ValueType>::initialize(const Vertex *const vertex) const {
  return std::make_shared<FunctionMemory<ValueType> >(EdgeDoFMacroVertexFunctionMemorySize, vertex->getNumNeighborEdges(), minLevel_,
                                                      maxLevel_);
}

template<typename ValueType>
std::shared_ptr<FunctionMemory<ValueType> >
EdgeDoFMacroEdgeFunctionMemoryDataHandling<ValueType>::initialize(const Edge *const edge) const {
  return std::make_shared<FunctionMemory<ValueType> >(EdgeDoFMacroEdgeFunctionMemorySize, edge->getNumNeighborFaces(), minLevel_,
                                                      maxLevel_);
}

template<typename ValueType>
std::shared_ptr<FunctionMemory<ValueType> > EdgeDoFMacroFaceFunctionMemoryDataHandling<ValueType>::initialize(const Face *const) const {
  return std::make_shared<FunctionMemory<ValueType> >(EdgeDoFMacroFaceFunctionMemorySize, 0, minLevel_, maxLevel_);
}

}/// namespace hhg