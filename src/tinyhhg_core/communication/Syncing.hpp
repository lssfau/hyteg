#pragma once

#include "core/DataTypes.h"

template < typename funcType >
class P2Function;

namespace hhg {
namespace communication {

using walberla::uint_t;

template < typename funcType >
void syncFunctionBetweenPrimitives( const funcType& function, const uint_t& level );

template < typename ValueType >
void syncP2FunctionBetweenPrimitives( const P2Function< ValueType >& function, const uint_t& level );

} // namespace communication
} // namespace hhg
