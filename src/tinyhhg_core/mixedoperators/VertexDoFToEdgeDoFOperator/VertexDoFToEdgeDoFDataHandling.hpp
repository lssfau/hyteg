#pragma once

#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/primitives/all.hpp"

namespace hhg {

class MacroEdgeVertexDoFToEdgeDoFDataHandling : public StencilMemoryDataHandling< StencilMemory < real_t >, Edge> {

public:

  MacroEdgeVertexDoFToEdgeDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel);
  inline std::shared_ptr< StencilMemory< real_t > > initialize(const Edge *const edge) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class MacroFaceVertexDoFToEdgeDoFDataHandling : public StencilMemoryDataHandling< StencilMemory < real_t >, Face> {

public:

  MacroFaceVertexDoFToEdgeDoFDataHandling(const uint_t &minLevel, const uint_t &maxLevel);
  inline std::shared_ptr< StencilMemory < real_t > > initialize(const Face *const face) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};


}//namespace hhg
