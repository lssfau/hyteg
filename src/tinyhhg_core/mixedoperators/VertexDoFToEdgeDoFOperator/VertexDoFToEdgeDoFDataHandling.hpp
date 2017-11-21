#pragma once

#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "VertexDoFToEdgeDoFMemory.hpp"

namespace hhg {

template<typename ValueType>
class EdgeVertexDoFToEdgeDoFStencilMemoryDataHandling
  : public StencilMemoryDataHandling<EdgeVertexDoFToEdgeDoFStencilMemory<ValueType>, Edge> {

public:

  EdgeVertexDoFToEdgeDoFStencilMemoryDataHandling(const uint_t &minLevel, const uint_t &maxLevel) : minLevel_(minLevel),
                                                                                                    maxLevel_(maxLevel) {}

  inline std::shared_ptr<EdgeVertexDoFToEdgeDoFStencilMemory<ValueType> > initialize(const Edge *const edge) const{
    return std::make_shared< EdgeVertexDoFToEdgeDoFStencilMemory< ValueType > >( EdgeVertexDoFToEdgeDoFStencilSize, edge->getNumNeighborFaces(), minLevel_, maxLevel_ );
  }

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template<typename ValueType>
class FaceVertexDoFToEdgeDoFStencilMemoryDataHandling
  : public StencilMemoryDataHandling<FaceVertexDoFToEdgeDoFStencilMemory<ValueType>, Face> {

public:

  FaceVertexDoFToEdgeDoFStencilMemoryDataHandling(const uint_t &minLevel, const uint_t &maxLevel) : minLevel_(minLevel),
                                                                                                    maxLevel_(maxLevel) {}

  inline std::shared_ptr<FaceVertexDoFToEdgeDoFStencilMemory<ValueType> > initialize(const Face *const face) const{
    return std::make_shared< FaceVertexDoFToEdgeDoFStencilMemory< ValueType > >( FaceVertexDoFToEdgeDoFStencilSize, 0, minLevel_, maxLevel_ );
  }

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};


}//namespace hhg
