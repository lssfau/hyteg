/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "Syncing.hpp"

#include "hyteg/functions/Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/egfunctionspace/EGFunction.hpp"
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

// Version for VectorFunctions
template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< vType >& vecFunc, const uint_t& level )
{
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
   {
      syncFunctionBetweenPrimitives( vecFunc[idx], level );
   }
}

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< vType >& vecFunc, const uint_t& level )
{
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
   {
      syncP2FunctionBetweenPrimitives( vecFunc[idx], level );
   }
}

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const dg::DGVectorFunction< vType >& vecFunc, const uint_t& level )
{
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
   {
      syncDGFunctionBetweenPrimitives( vecFunc[idx], level );
   }
}

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const EGFunction< vType >& p1dgeFunc, const uint_t& level )
{
   auto& vecFunc = *(p1dgeFunc.getConformingPart());
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
      syncFunctionBetweenPrimitives< hyteg::vertexdof::VertexDoFFunction< vType > >( vecFunc[idx], level );
}

template void syncP2FunctionBetweenPrimitives( const P2Function< double >& function, const uint_t& level );
template void syncP2FunctionBetweenPrimitives( const P2Function< float >& function, const uint_t& level );
template void syncP2FunctionBetweenPrimitives( const P2Function< int32_t >& function, const uint_t& level );
template void syncP2FunctionBetweenPrimitives( const P2Function< int64_t >& function, const uint_t& level );

template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< double >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< float >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< int32_t >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< int64_t >& function, const uint_t& level );

template void syncFunctionBetweenPrimitives( const EdgeDoFFunction< double >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction< float >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction< int32_t >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction< int64_t >& function, const uint_t& level );

template void syncFunctionBetweenPrimitives( const P2Function< double >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const P2Function< float >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const P2Function< int32_t >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const P2Function< int64_t >& function, const uint_t& level );

template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< double >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< float >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< int32_t >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< int64_t >& function, const uint_t& level );

template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< double >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< float >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< int32_t >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< int64_t >& function, const uint_t& level );

template void syncVectorFunctionBetweenPrimitives( const EGFunction< double >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const EGFunction< int32_t >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const EGFunction< int64_t >& function, const uint_t& level );

} // namespace communication
} // namespace hyteg
