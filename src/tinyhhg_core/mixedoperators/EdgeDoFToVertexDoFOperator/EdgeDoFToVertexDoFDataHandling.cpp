#include "EdgeDoFToVertexDoFDataHandling.hpp"

namespace hhg {


////////// Macro Vertex DataHandling ///////////
MacroVertexEdgeDoFToVertexDoFDataHandling::MacroVertexEdgeDoFToVertexDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel)
  :minLevel_(minLevel),
   maxLevel_(maxLevel)
{}

inline std::shared_ptr<StencilMemory<real_t> > MacroVertexEdgeDoFToVertexDoFDataHandling::initialize(const Vertex *const vertex) const {
  return std::make_shared< StencilMemory< real_t > >( macroVertexVertexDoFToEdgeDoFStencilSize,
                                                      vertex->getNumNeighborEdges(),
                                                      minLevel_,
                                                      maxLevel_ );
}

////////// Macro Edge DataHandling //////////
MacroEdgeEdgeDoFToVertexDoFDataHandling::MacroEdgeEdgeDoFToVertexDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel)
  :minLevel_(minLevel),
   maxLevel_(maxLevel)
{}

inline std::shared_ptr<StencilMemory<real_t> > MacroEdgeEdgeDoFToVertexDoFDataHandling::initialize(const Edge *const edge) const {
  return std::make_shared< StencilMemory< real_t > >( macroEdgeVertexDoFToEdgeDoFStencilSize,
                                                      edge->getNumNeighborFaces(),
                                                      minLevel_,
                                                      maxLevel_ );
}

////////// Macro Face DataHandling //////////
MacroFaceEdgeDoFToVertexDoFDataHandling::MacroFaceEdgeDoFToVertexDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel)
  :minLevel_(minLevel),
   maxLevel_(maxLevel)
{}

inline std::shared_ptr<StencilMemory<real_t> > MacroFaceEdgeDoFToVertexDoFDataHandling::initialize(const Face *const face) const {
  WALBERLA_UNUSED(face);
  return std::make_shared< StencilMemory< real_t > >( macroFaceVertexDoFToEdgeDoFStencilSize,
                                                      0,
                                                      minLevel_,
                                                      maxLevel_ );
}


////////// Stencil sizes //////////
uint_t macroVertexVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 2 * numDependencies;
}

uint_t macroEdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  return 2 + 5 * numDependencies;
}

uint_t macroFaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies) {
  WALBERLA_UNUSED(level);
  WALBERLA_UNUSED(numDependencies);
  return 12;
}

} //namespace hhg

