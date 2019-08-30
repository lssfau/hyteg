#include "Syncing.hpp"

#include "hyteg/Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {
namespace communication {

using walberla::uint_t;

template < typename funcType >
void syncFunctionBetweenPrimitives( const funcType& function, const uint_t& level )
{
   function.template communicate< Vertex, Edge >( level );
   function.template communicate< Edge, Face >( level );
   function.template communicate< Face, Cell >( level );

   function.template communicate< Cell, Face >( level );
   function.template communicate< Face, Edge >( level );
   function.template communicate< Edge, Vertex >( level );
}

template < typename ValueType >
void syncP2FunctionBetweenPrimitives( const P2Function< ValueType >& function, const uint_t& level )
{
   syncFunctionBetweenPrimitives< hyteg::vertexdof::VertexDoFFunction< ValueType > >( function.getVertexDoFFunction(), level );
   syncFunctionBetweenPrimitives< hyteg::EdgeDoFFunction< ValueType > >( function.getEdgeDoFFunction(), level );
}

template void syncP2FunctionBetweenPrimitives( const P2Function< double >& function, const uint_t& level );

template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction<double>& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction<int>& function, const uint_t& level );

template void syncFunctionBetweenPrimitives( const EdgeDoFFunction<double>& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction<int>& function, const uint_t& level );

//template void syncFunctionBetweenPrimitives( const P1Function<double>& function, const uint_t& level );

} // namespace communication
} // namespace hyteg