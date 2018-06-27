#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitives/all.hpp"

namespace hhg {
namespace communication {

using walberla::uint_t;

template < typename funcType >
void syncFunctionBetweenPrimitives( const funcType& function, const uint_t level )
{
   function.template communicate< Vertex, Edge >( level );
   function.template communicate< Edge, Face >( level );
   function.template communicate< Face, Cell >( level );

   function.template communicate< Cell, Face >( level );
   function.template communicate< Face, Edge >( level );
   function.template communicate< Edge, Vertex >( level );
}

template < typename ValueType >
void syncP2FunctionBetweenPrimitives( const P2Function< ValueType >& function, const uint_t level )
{
   syncFunctionBetweenPrimitives< hhg::vertexdof::VertexDoFFunction< ValueType > >( *function.getVertexDoFFunction(), level );
   syncFunctionBetweenPrimitives< hhg::EdgeDoFFunction< ValueType > >( *function.getEdgeDoFFunction(), level );
}

} // namespace communication
} // namespace hhg
