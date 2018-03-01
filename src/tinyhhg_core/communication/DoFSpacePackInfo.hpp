#pragma once

#include "PackInfo.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"

//#include "tinyhhg_core/primitives/all.hpp"


namespace hhg{

namespace communication{

template< typename ValueType >
class DoFSpacePackInfo : public communication::PackInfo
{
public:

  DoFSpacePackInfo(uint_t level,
                   PrimitiveDataID <FunctionMemory<ValueType>, Vertex> dataIDVertex,
                   PrimitiveDataID <FunctionMemory<ValueType>, Edge> dataIDEdge,
                   PrimitiveDataID <FunctionMemory<ValueType>, Face> dataIDFace,
                   std::weak_ptr <PrimitiveStorage> storage)
    : level_(level),
      dataIDVertex_(dataIDVertex),
      dataIDEdge_(dataIDEdge),
      dataIDFace_(dataIDFace),
      dataIDCell_(storage.lock()->generateInvalidPrimitiveDataID< FunctionMemory< ValueType >, Cell >()),
      storage_(storage) {}

  DoFSpacePackInfo(uint_t level,
                   PrimitiveDataID <FunctionMemory<ValueType>, Vertex> dataIDVertex,
                   PrimitiveDataID <FunctionMemory<ValueType>, Edge> dataIDEdge,
                   PrimitiveDataID <FunctionMemory<ValueType>, Face> dataIDFace,
                   PrimitiveDataID <FunctionMemory<ValueType>, Cell> dataIDCell,
                   std::weak_ptr <PrimitiveStorage> storage)
    : level_(level),
      dataIDVertex_(dataIDVertex),
      dataIDEdge_(dataIDEdge),
      dataIDFace_(dataIDFace),
      dataIDCell_(dataIDCell),
      storage_(storage) {}

protected:

  uint_t level_;
  PrimitiveDataID<FunctionMemory< ValueType >, Vertex> dataIDVertex_;
  PrimitiveDataID<FunctionMemory< ValueType >, Edge> dataIDEdge_;
  PrimitiveDataID<FunctionMemory< ValueType >, Face> dataIDFace_;
  PrimitiveDataID<FunctionMemory< ValueType >, Cell> dataIDCell_;
  std::weak_ptr<hhg::PrimitiveStorage> storage_;
};


}// namespace communication
}// namesapce hhg
