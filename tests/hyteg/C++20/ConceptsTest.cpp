/*
* Copyright (c) 2024 Marcus Mohr.
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

// check pipeline compilers for support of C++20 concepts

#include <concepts>
#include <iostream>

#include "core/Environment.h"

#include "hyteg/functions/CSFVectorFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/types/Concepts.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

// Test template
template< typename T >
requires std::integral<T>
struct myClass {
  T attribute;
};

template< concepts::concrete_primitive primitive_t >
void accessTypeInfo( const primitive_t& primitive ) {
  WALBERLA_LOG_INFO_ON_ROOT( "" << "Type is " << primitive.getType() );
}

// Test using concept
void getFunctionName( const concepts::fe_function auto& function ) {
  std::cout << "Name of function is ' " << function.getFunctionName() << "'" << std::endl;
}

int main( int argc, char** argv )
{

   // environment stuff
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   myClass< uint > obj1;

   obj1.attribute = 2UL;

   WALBERLA_UNUSED( obj1 );

   MeshInfo                            mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P2Function< real_t > p2Func( "Lagrange, 2nd order", storage, 2, 2 );
   P1VectorFunction< real_t > p1Func( "P1 Vector Function", storage, 2, 2 );

   getFunctionName( p1Func );
   getFunctionName( p2Func );
   // getFunctionName( obj1 );

   WALBERLA_LOG_INFO_ON_ROOT( "If you see this message, compiling with concepts worked ;-)" );

   return EXIT_SUCCESS;
}
