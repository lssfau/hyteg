#include "Syncing.hpp"

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/primitives/all.hpp"

namespace hhg {
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
   syncFunctionBetweenPrimitives< hhg::vertexdof::VertexDoFFunction< ValueType > >( function.getVertexDoFFunction(), level );
   syncFunctionBetweenPrimitives< hhg::EdgeDoFFunction< ValueType > >( function.getEdgeDoFFunction(), level );
}

template void syncP2FunctionBetweenPrimitives( const P2Function< float >& function, const uint_t& level );
template void syncP2FunctionBetweenPrimitives( const P2Function< double >& function, const uint_t& level );

template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction<float>& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction<double>& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction<int>& function, const uint_t& level );

template void syncFunctionBetweenPrimitives( const EdgeDoFFunction<float>& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction<double>& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction<int>& function, const uint_t& level );

//template void syncFunctionBetweenPrimitives( const P1Function<float>& function, const uint_t& level );
//template void syncFunctionBetweenPrimitives( const P1Function<double>& function, const uint_t& level );

} // namespace communication
} // namespace hhg