/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/forms/form_fenics_generated/p1_tet_diffusion.h"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "constantStencilOperator/P1ConstantOperator.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

/// This is rather an instructive example illustrating a cool property of std::shared_ptr< void >.
/// When std::shared_ptr< void > is constructed from another pointer of type T*, the type is erased but the deleter is stored!
/// This means although we perform type erasure, the allocated memory is safely deleted.
/// The exact same mechanism is used in HyTeG's PrimitiveData class to attach arbitrary simulation data to Primitives.
class DeleterTest
{
 public:
   DeleterTest() { WALBERLA_LOG_INFO_ON_ROOT( "Constructor of \"DeleterTest\" is called." ) }

   ~DeleterTest() { WALBERLA_LOG_INFO_ON_ROOT( "Destructor of \"DeleterTest\" is called." ) }

   void sayHello() { WALBERLA_LOG_INFO_ON_ROOT( "Hello!" ) }
};

// Pure C++.
void testSharedPtrDeleter()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Pure C++ example." )

   std::shared_ptr< void > voidPtr( new DeleterTest() );

   static_cast< DeleterTest* >( voidPtr.get() )->sayHello();
}

// Using PrimitiveData.
void testPrimitiveData()
{
   WALBERLA_LOG_INFO_ON_ROOT( "HyTeG example." )

   internal::PrimitiveData pdata( std::make_shared< DeleterTest >() );

   pdata.get< DeleterTest >()->sayHello();
}

/// --- Actual tests --------------------------------------------------------------------------------------------------------

/// This function tests the destructor, copy-constructor and copy-assignment for various functions.
/// The issue is that functions do not "own" the memory they allocate.
/// They only carry handles (PrimitiveDataIDs) that can be used to retrieve pointers from Primitives.
/// Having the Primitives own the actual allocated memory allows for "straightforward" data migration in parallel (think load balancing).
/// Without any information about functions and operators, the PrimitiveStorage can migrate the data to other processes.
/// On the downside, copying or deleting functions means that the memory that is attached to the Primitives must be handled accordingly.
///
/// Eventually functions are only handles to the memory. Copying does not copy the allocate memory.
/// However, we want the memory to be deallocated once all handles are destroyed. This is tested here.
///
/// See issue https://i10git.cs.fau.de/hyteg/hyteg/-/issues/186
void testRuleOfThree()
{
   uint_t level = 2;

   auto                  meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_24el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto            storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto primitiveIDs = storage->getPrimitiveIDs();

   ///////////////
   // VertexDoF //
   ///////////////

   for ( auto id : primitiveIDs )
   {
      auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
      WALBERLA_CHECK_EQUAL( numDataEntries, 0, "No data should be allocated before creating functions." );
   }

   {
      vertexdof::VertexDoFFunction< real_t > f( "f", storage, level, level );

      auto f_copy( f );
      auto f_copy_2 = f;
      WALBERLA_CHECK_EQUAL( f.getFunctionName(), f_copy_2.getFunctionName() );

      f             = f_copy_2;

      WALBERLA_LOG_INFO_ON_ROOT( "Just printing to avoid stuff getting optimized away... " << f.getFunctionName() << ", "
                                                                                           << f_copy.getFunctionName() << ", "
                                                                                           << f_copy_2.getFunctionName() );

      for ( auto id : primitiveIDs )
      {
         auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
         WALBERLA_CHECK_EQUAL( numDataEntries, 1, "Data should be allocated after creating VertexDoFFunction." );
      }
   }

   for ( auto id : primitiveIDs )
   {
      auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
      WALBERLA_CHECK_EQUAL( numDataEntries, 0, "No data should be allocated after leaving scope." );
   }

   /////////////
   // EdgeDoF //
   /////////////

   for ( auto id : primitiveIDs )
   {
      auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
      WALBERLA_CHECK_EQUAL( numDataEntries, 0, "No data should be allocated before creating functions." );
   }

   {
      EdgeDoFFunction< real_t > f( "f", storage, level, level );

      auto f_copy( f );
      auto f_copy_2 = f;
      f             = f_copy_2;

      WALBERLA_LOG_INFO_ON_ROOT( "Just printing to avoid stuff getting optimized away... " << f.getFunctionName() << ", "
                                                                                           << f_copy.getFunctionName() << ", "
                                                                                           << f_copy_2.getFunctionName() );

      for ( auto id : primitiveIDs )
      {
         auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
         WALBERLA_CHECK_EQUAL( numDataEntries, 1, "Data should be allocated after creating EdgeDoFFunction." );
      }
   }

   for ( auto id : primitiveIDs )
   {
      auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
      WALBERLA_CHECK_EQUAL( numDataEntries, 0, "No data should be allocated after leaving scope." );
   }

   ///////////////
   // VolumeDoF //
   ///////////////

   for ( auto id : primitiveIDs )
   {
      auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
      WALBERLA_CHECK_EQUAL( numDataEntries, 0, "No data should be allocated before creating functions." );
   }

   {
      volumedofspace::VolumeDoFFunction< real_t > f(
          "f", storage, level, level, 2, volumedofspace::indexing::VolumeDoFMemoryLayout::AoS );

      auto f_copy( f );
      auto f_copy_2 = f;
      f             = f_copy_2;

      WALBERLA_LOG_INFO_ON_ROOT( "Just printing to avoid stuff getting optimized away... " << f.getFunctionName() << ", "
                                                                                           << f_copy.getFunctionName() << ", "
                                                                                           << f_copy_2.getFunctionName() );

      for ( auto id : primitiveIDs )
      {
         if ( storage->cellExistsLocally( id ) )
         {
            auto numDataEntries = storage->getCell( id )->getNumberOfDataEntries();
            WALBERLA_CHECK_EQUAL( numDataEntries, 5, "Data should be allocated after creating VolumeDoFFunction." );
         }
         else
         {
            auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
            WALBERLA_CHECK_EQUAL( numDataEntries, 0, "Data should be allocated after creating VolumeDoFFunction." );
         }
      }
   }

   for ( auto id : primitiveIDs )
   {
      auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
      WALBERLA_CHECK_EQUAL( numDataEntries, 0, "No data should be allocated after leaving scope." );
   }

   //////////////////
   /// Composites ///
   //////////////////

   for ( auto id : primitiveIDs )
   {
      auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
      WALBERLA_CHECK_EQUAL( numDataEntries, 0, "No data should be allocated before creating functions." );
   }

   {
      auto                     basis = std::make_shared< dg::DGBasisLinearLagrange_Example >();
      dg::DGFunction< real_t > f_dg( "f_dg", storage, level, level, basis, 1 );
      WALBERLA_LOG_INFO_ON_ROOT( "Just printing to avoid stuff getting optimized away... " << f_dg.getFunctionName() );

      P2Function< real_t > f_p2( "f_p2", storage, level, level );
      WALBERLA_LOG_INFO_ON_ROOT( "Just printing to avoid stuff getting optimized away... " << f_p2.getFunctionName() );

      P2P1TaylorHoodFunction< real_t > f_p2p1( "f_p2p1", storage, level, level );
      WALBERLA_LOG_INFO_ON_ROOT( "Just printing to avoid stuff getting optimized away... " << f_p2p1.getFunctionName() );

      for ( auto id : primitiveIDs )
      {
         auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
         WALBERLA_CHECK_GREATER( numDataEntries, 0, "Data should be allocated after creating VolumeDoFFunction." );
      }
   }

   for ( auto id : primitiveIDs )
   {
      auto numDataEntries = storage->getPrimitive( id )->getNumberOfDataEntries();
      WALBERLA_CHECK_EQUAL( numDataEntries, 0, "No data should be allocated after leaving scope." );
   }
}

template< typename func_t >
void testCopyAssignment() {

   WALBERLA_LOG_INFO_ON_ROOT( "**** " << FunctionTrait< func_t >::getTypeName()  << " ****" );

   uint_t level = 2;

   auto                  meshInfo1 = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage1( meshInfo1, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto            storage1 = std::make_shared< PrimitiveStorage >( setupStorage1, 1 );
 
   func_t orig( "original", storage1, level, level );

   auto                  meshInfo2 = MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" );
   SetupPrimitiveStorage setupStorage2( meshInfo2, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto            storage2 = std::make_shared< PrimitiveStorage >( setupStorage2, 1 );

   func_t clone( "should become 'original'", storage2, level, level );

   clone = orig;
   WALBERLA_CHECK_EQUAL( clone.getStorage(), orig.getStorage() );
   WALBERLA_CHECK_EQUAL( clone.getFunctionName(), orig.getFunctionName() );
   WALBERLA_CHECK_EQUAL( clone.getMinLevel(), orig.getMinLevel() );
   WALBERLA_CHECK_EQUAL( clone.getMaxLevel(), orig.getMaxLevel() );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------------------------" )
   WALBERLA_LOG_INFO_ON_ROOT( "--- std::shared_ptr< void > deleter examples ---" )
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------------------------" )
   testSharedPtrDeleter();
   testPrimitiveData();
   WALBERLA_LOG_INFO_ON_ROOT( "------------------------------------------------" )
   testRuleOfThree();

   WALBERLA_LOG_INFO_ON_ROOT( "===============================" );
   WALBERLA_LOG_INFO_ON_ROOT( "==> Testing Copy Assignment <==" );
   WALBERLA_LOG_INFO_ON_ROOT( "===============================" );

   testCopyAssignment< P1Function< real_t > >();
   testCopyAssignment< P2Function< real_t > >();
   testCopyAssignment< EdgeDoFFunction< real_t > >();

   testCopyAssignment< P1VectorFunction< real_t > >();
   testCopyAssignment< P2VectorFunction< real_t > >();

   testCopyAssignment< P0Function< real_t > >();
   testCopyAssignment< P2P1TaylorHoodFunction< real_t > >();

   return 0;
}
