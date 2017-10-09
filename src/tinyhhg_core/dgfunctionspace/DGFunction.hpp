#pragma once

#include "tinyhhg_core/Function.hpp"
#include "DGPackInfo.hpp"

using namespace hhg;

template< typename ValueType >
class DGFunction : public Function< DGFunction< ValueType > > {
public:

  using Function<DGFunction<ValueType> >::storage_;
  using Function<DGFunction<ValueType> >::communicators_;

  DGFunction(const std::string &name, const std::shared_ptr<PrimitiveStorage> &storage, uint_t minLevel,
             uint_t maxLevel) :
      Function<DGFunction<ValueType> >(name, storage, minLevel, maxLevel)
  {
    auto faceDGFunctionMemoryDataHandling =
        std::make_shared<FaceDGFunctionMemoryDataHandling < ValueType> > (minLevel, maxLevel);
    auto edgeDGFunctionMemoryDataHandling =
        std::make_shared<EdgeDGFunctionMemoryDataHandling < ValueType> > (minLevel, maxLevel);
    auto vertexDGFunctionMemoryDataHandling =
        std::make_shared<VertexDGFunctionMemoryDataHandling < ValueType> > (minLevel, maxLevel);
    storage->addFaceData(faceDataID_, faceDGFunctionMemoryDataHandling, name);
    storage->addEdgeData(edgeDataID_, edgeDGFunctionMemoryDataHandling, name);
    storage->addVertexData(vertexDataID_, vertexDGFunctionMemoryDataHandling, name);
    for (uint_t level = minLevel; level <= maxLevel; ++level)
    {
      communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
      communicators_[level]->addPackInfo(
          std::make_shared<DGPackInfo<ValueType> >(level, vertexDataID_, edgeDataID_, faceDataID_, storage_));
    }
  }

  const PrimitiveDataID<FunctionMemory<ValueType>, Vertex> &getVertexDataID() const
  { return vertexDataID_; }

  const PrimitiveDataID<FunctionMemory<ValueType>, Edge> &getEdgeDataID() const
  { return edgeDataID_; }

  const PrimitiveDataID<FunctionMemory<ValueType>, Face> &getFaceDataID() const
  { return faceDataID_; }

};
