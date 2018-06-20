#pragma once

#include <string>

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"


namespace hhg {
namespace communication {

using walberla::uint_t;


template < typename funcType >
void syncFunctionBetweenPrimitives( const funcType& function, uint_t level )
{
   function.template communicate< Vertex, Edge >( level );
   function.template communicate< Edge, Face >( level );
   function.template communicate< Face, Cell >( level );

   function.template communicate< Cell, Face >( level );
   function.template communicate< Face, Edge >( level );
   function.template communicate< Edge, Vertex >( level );
}

template < typename ValueType >
void syncFunctionBetweenPrimitives( P2Function< ValueType > function, uint_t level)
{
   syncFunctionBetweenPrimitives< hhg::vertexdof::VertexDoFFunction< ValueType > >(*function.getVertexDoFFunction(), level);
   syncFunctionBetweenPrimitives< hhg::EdgeDoFFunction< ValueType > >(*function.getEdgeDoFFunction(), level);
}

} // namespace communication
} // namespace hhg
