#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/communication/PackInfo.hpp"

namespace hyteg {

class Vertex;
class Edge;
class Face;
class Cell;
class PrimitiveStorage;

template < typename ValueType >
class FunctionMemory;
template < typename DataType, typename PrimitiveType >
class PrimitiveDataID;

namespace communication {

using walberla::uint_t;

template < typename ValueType >
class DoFSpacePackInfo : public communication::PackInfo
{
 public:
   DoFSpacePackInfo( uint_t                                                 level,
                     PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                     PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                     PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                     std::weak_ptr< PrimitiveStorage >                      storage );

   DoFSpacePackInfo( uint_t                                                 level,
                     PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                     PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                     PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                     PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                     std::weak_ptr< PrimitiveStorage >                      storage );

 protected:
   uint_t                                                 level_;
   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell_;
   std::weak_ptr< hyteg::PrimitiveStorage >                 storage_;
};

} // namespace communication
} // namespace hyteg
