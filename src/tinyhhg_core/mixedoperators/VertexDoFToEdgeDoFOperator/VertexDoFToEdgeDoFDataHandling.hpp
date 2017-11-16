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

  inline std::shared_ptr<EdgeVertexDoFToEdgeDoFStencilMemory<ValueType> > initialize(const Edge *const edge) const;

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

  inline std::shared_ptr<FaceVertexDoFToEdgeDoFStencilMemory<ValueType> > initialize(const Face *const face) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
std::shared_ptr< EdgeVertexDoFToEdgeDoFStencilMemory< ValueType > > EdgeVertexDoFToEdgeDoFStencilMemoryDataHandling< ValueType >::initialize( const Edge * const edge ) const
{
  return std::make_shared< EdgeVertexDoFToEdgeDoFStencilMemory< ValueType > >( EdgeVertexDoFToEdgeDoFStencil, edge->getNumNeighborFaces(), minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< FaceVertexDoFToEdgeDoFStencilMemory< ValueType > > FaceVertexDoFToEdgeDoFStencilMemoryDataHandling< ValueType >::initialize( const Face * const face ) const
{
  return std::make_shared< FaceVertexDoFToEdgeDoFStencilMemory< ValueType > >( FaceVertexDoFToEdgeDoFStencil, 0, minLevel_, maxLevel_ );
}

}//namespace hhg
