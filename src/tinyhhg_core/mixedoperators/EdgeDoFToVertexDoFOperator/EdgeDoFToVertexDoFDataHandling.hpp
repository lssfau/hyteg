#pragma once

#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/primitives/all.hpp"

namespace hhg {


class MacroVertexEdgeDoFToVertexDoFDataHandling : public StencilMemoryDataHandling<StencilMemory < real_t>, Vertex>{

public:
  MacroVertexEdgeDoFToVertexDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel);
  std::shared_ptr<StencilMemory< real_t > > initialize(const Vertex *const vertex) const final;

private:
  uint_t minLevel_;
  uint_t maxLevel_;
};

class MacroEdgeEdgeDoFToVertexDoFDataHandling : public StencilMemoryDataHandling<StencilMemory < real_t>, Edge> {

public:
  MacroEdgeEdgeDoFToVertexDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel);
  std::shared_ptr<StencilMemory< real_t > > initialize(const Edge *const edge) const final;

private:
  uint_t minLevel_;
  uint_t maxLevel_;
};

class MacroFaceEdgeDoFToVertexDoFDataHandling : public StencilMemoryDataHandling<StencilMemory < real_t>, Face> {

public:
  MacroFaceEdgeDoFToVertexDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel);
  std::shared_ptr<StencilMemory< real_t > > initialize(const Face *const face) const final;

private:
  uint_t minLevel_;
  uint_t maxLevel_;
};

/// \param level the stencil size is independent of the level
/// \param numDependencies number of adjacent edges of the vertex
/// \return number of entries in the stencil on a macro vertex, be aware that this one entry to large on the boundaries
inline uint_t macroVertexVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies);

/// \param level the stencil size is independent of the level
/// \param numDependencies number of adjacent faces of the edge
/// \return number of entries in the stencil on a macro edge
inline uint_t macroEdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies);

/// \param level the stencil size is independent of the level
/// \param numDependencies on a macro face this is always constant
/// \return number of entries in the stencil on a macro face
inline uint_t macroFaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const uint_t &numDependencies);

}//namespace hhg