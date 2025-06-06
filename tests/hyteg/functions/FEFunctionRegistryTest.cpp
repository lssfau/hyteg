/*
 * Copyright (c) 2023-2025 Marcus Mohr.
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
#include "hyteg/functions/FEFunctionRegistry.hpp"

#include <cstdint>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"

#include "hyteg/composites/CCRStokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// using walberla::uint_t;
using namespace hyteg;
using namespace hyteg::functionTraits;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t minLevel{ 1 };
   const uint_t maxLevel{ 2 };

   // Prepare mesh and primitive storages
   MeshInfo                            mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "3D/pyramid_2el.msh" ) );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage   = std::make_shared< PrimitiveStorage >( setupStorage );
   std::shared_ptr< PrimitiveStorage > storageDG = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // Setup some functions
   P0Function< real_t >  p0ScalarFunc1( "P0 scalar function 1", storageDG, minLevel, maxLevel );
   P0Function< int64_t > p0ScalarFunc2( "P0 scalar function 2", storageDG, minLevel, maxLevel );

   P1Function< real_t >  p1ScalarFunc1( "P1 scalar function 1", storage, minLevel, maxLevel );
   P1Function< int32_t > p1ScalarFunc2( "P1 scalar function 2", storage, minLevel, maxLevel );

   P2Function< real_t >  p2ScalarFunc1( "P2 scalar function 1", storage, minLevel, maxLevel );
   P2Function< int64_t > p2ScalarFunc2( "P2 scalar function 2", storage, minLevel, maxLevel );

   DG1Function< real_t >  dg1ScalarFunc1( "DG1 scalar function 1", storageDG, minLevel, maxLevel );
   DG1Function< int32_t > dg1ScalarFunc2( "DG1 scalar function 2", storageDG, minLevel, maxLevel );

   P1VectorFunction< real_t >  p1VectorFunc( "P1 vector function", storage, minLevel, maxLevel );
   P2VectorFunction< int64_t > p2VectorFunc( "P2 vector function", storage, minLevel, maxLevel );

   P2P1TaylorHoodFunction< real_t > stokesFunc( "Stokes function", storage, minLevel, maxLevel );
   CCRStokesFunction< real_t >      ccrStokesFunc( "CCR Stokes function", storage, minLevel, maxLevel );

   EGFunction< real_t > egVectorFunc( "EG vector function", storageDG, minLevel, maxLevel );

   n1e1::N1E1VectorFunction< real_t > n1e1VectorFunc( "N1E1 vector function", storage, minLevel, maxLevel );

   P2PlusBubbleFunction< real_t > p2BubbleFunc( "bubble enhance P2 function", storageDG, minLevel, maxLevel );

   // create a registry and register some functions
   FEFunctionRegistry registry;

   registry.add( p0ScalarFunc1 );
   registry.add( p0ScalarFunc2 );

   registry.add( p1ScalarFunc1 );
   registry.add( p1ScalarFunc2 );

   registry.add( p2ScalarFunc1 );

   // check name extraction
   std::vector< std::string > names;
   registry.extractFunctionNames( names, functionTraits::P2_FUNCTION );
   WALBERLA_CHECK_EQUAL( names[0], "P2 scalar function 1" );

   registry.extractFunctionNames( names, functionTraits::DG_FUNCTION );
   WALBERLA_CHECK_EQUAL( names[1], "P0 scalar function 1" );
   WALBERLA_CHECK_EQUAL( names[2], "P0 scalar function 2" );

   // add further functions
   registry.add( dg1ScalarFunc1 );
   registry.add( egVectorFunc );
   registry.add( n1e1VectorFunc );
   registry.add( p1VectorFunc );
   registry.add( p2VectorFunc );
   registry.add( p2BubbleFunc );

   // adding this will increate the count of P1Function and P2VectorFunction
   registry.add( stokesFunc );

   // adding this will increate the count of P2PlusBubbleVectorFunction and DG1Function
   registry.add( ccrStokesFunc );

   // check name extraction individually for all kinds
   std::map< functionTraits::FunctionKind, uint_t > kindCount{ { DG_FUNCTION, 4 },
                                                               { P1_FUNCTION, 3 },
                                                               { P2_FUNCTION, 1 },
                                                               { P2_PLUS_BUBBLE_FUNCTION, 1 },
                                                               { P2_PLUS_BUBBLE_VECTOR_FUNCTION, 1 },
                                                               { P1_VECTOR_FUNCTION, 1 },
                                                               { P2_VECTOR_FUNCTION, 2 },
                                                               { N1E1_VECTOR_FUNCTION, 1 },
                                                               { EG_FUNCTION, 1 } };

   for ( const auto& entry : kindCount )
   {
      names.clear();
      registry.extractFunctionNames( names, entry.first );
      WALBERLA_CHECK_EQUAL( names.size(), entry.second );
   }

   // selectively check function extraction
   const auto p1Funcs = registry.getP1Functions();
   WALBERLA_CHECK_EQUAL( p1Funcs.size(), 3 );

   const auto p2VecFuncs = registry.getP2VectorFunctions();
   WALBERLA_CHECK_EQUAL( p2VecFuncs.size(), 2 );

   const auto nedelecFuncs = registry.getN1E1VectorFunctions();
   WALBERLA_CHECK_EQUAL( nedelecFuncs.size(), 1 );

   // now go the full way to the function
   const auto p2Funcs = registry.getP2Functions();
   WALBERLA_CHECK_EQUAL( p2Funcs.size(), 1 );
   std::vector< P2Function< real_t > > realP2Funcs = p2Funcs.getFunctions< real_t >();
   WALBERLA_CHECK_EQUAL( p2ScalarFunc1.getFunctionName(), realP2Funcs[0].getFunctionName() );

   // test deregistering a function works
   registry.remove( p2VectorFunc );
   WALBERLA_CHECK_EQUAL( registry.getP2VectorFunctions().size(), --kindCount[P2_VECTOR_FUNCTION] );

   registry.remove( p0ScalarFunc1 );
   WALBERLA_CHECK_EQUAL( registry.getDGFunctions().size(), --kindCount[DG_FUNCTION] );

   registry.remove( stokesFunc );
   WALBERLA_CHECK_EQUAL( registry.getP2VectorFunctions().size(), --kindCount[P2_VECTOR_FUNCTION] );
   WALBERLA_CHECK_EQUAL( registry.getP1Functions().size(), --kindCount[P1_FUNCTION] );

   registry.remove( ccrStokesFunc );
   WALBERLA_CHECK_EQUAL( registry.getP2PlusBubbleVectorFunctions().size(), --kindCount[P2_PLUS_BUBBLE_VECTOR_FUNCTION] );
   WALBERLA_CHECK_EQUAL( registry.getDGFunctions().size(), --kindCount[DG_FUNCTION] );

   return EXIT_SUCCESS;
}
