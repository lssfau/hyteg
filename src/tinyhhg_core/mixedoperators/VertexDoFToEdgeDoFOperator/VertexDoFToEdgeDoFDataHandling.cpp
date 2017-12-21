#include "VertexDoFToEdgeDoFDataHandling.hpp"

namespace hhg{

///
/// \param level stencil size is independent of level
/// \param numDependencies number of adjacent faces of the edge
/// \return number of the stencil entries
inline uint_t macroEdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies)
{
  WALBERLA_UNUSED( level );
  return 2 + numDependencies;
}

inline uint_t macroFaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies)
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( numDependencies );
  return 4;
}

MacroEdgeVertexDoFToEdgeDoFDataHandling::MacroEdgeVertexDoFToEdgeDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel)
  :minLevel_(minLevel),
   maxLevel_(maxLevel)
{}

std::shared_ptr<StencilMemory< real_t >> MacroEdgeVertexDoFToEdgeDoFDataHandling::initialize(const Edge *const edge) const {
  return std::make_shared<StencilMemory<real_t> >(macroEdgeVertexDoFToEdgeDoFStencilSize,
                                                  edge->getNumNeighborFaces(),
                                                  minLevel_,
                                                  maxLevel_);
}

MacroFaceVertexDoFToEdgeDoFDataHandling::MacroFaceVertexDoFToEdgeDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel)
  :minLevel_(minLevel),
   maxLevel_(maxLevel)
{}

std::shared_ptr<StencilMemory< real_t >> MacroFaceVertexDoFToEdgeDoFDataHandling::initialize(const Face *const face) const {
  return std::make_shared<StencilMemory< real_t > >(macroFaceVertexDoFToEdgeDoFStencilSize,
                                                  face->getNumNeighborFaces(),
                                                  minLevel_,
                                                  maxLevel_);
}


} /// namespace hhg