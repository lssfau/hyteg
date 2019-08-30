#include "DoFSpacePackInfo.hpp"

#include "hyteg/FunctionMemory.hpp"
#include "hyteg/communication/PackInfo.hpp"
#include "hyteg/primitivedata/PrimitiveDataID.hpp"

namespace hyteg {
namespace communication {

template < typename ValueType >
DoFSpacePackInfo< ValueType >::DoFSpacePackInfo( uint_t                                                 level,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                                                 std::weak_ptr< PrimitiveStorage >                      storage )
: level_( level )
, dataIDVertex_( dataIDVertex )
, dataIDEdge_( dataIDEdge )
, dataIDFace_( dataIDFace )
, dataIDCell_( storage.lock()->generateInvalidPrimitiveDataID< FunctionMemory< ValueType >, Cell >() )
, storage_( storage )
{}

template < typename ValueType >
DoFSpacePackInfo< ValueType >::DoFSpacePackInfo( uint_t                                                 level,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                                                 PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                                                 std::weak_ptr< PrimitiveStorage >                      storage )
: level_( level )
, dataIDVertex_( dataIDVertex )
, dataIDEdge_( dataIDEdge )
, dataIDFace_( dataIDFace )
, dataIDCell_( dataIDCell )
, storage_( storage )
{}

template class DoFSpacePackInfo< double >;
template class DoFSpacePackInfo< int >;
template class DoFSpacePackInfo< uint_t >;

} // namespace communication
} // namespace hyteg