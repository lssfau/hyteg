#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"

#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"

using walberla::uint_t;

namespace hyteg {

template < typename FunctionTag_T >
void testFunctionProperties( const std::string& meshFileName,
                             const uint_t&      level,
                             const uint_t&      expectedGlobalDoFs,
                             const uint_t&      expectedInnerDoFs )
{
   auto storage = PrimitiveStorage::createFromGmshFile( meshFileName );

   auto storageInfo = storage->getGlobalInfo();

   const uint_t numGlobalDoFs = numberOfGlobalDoFs< FunctionTag_T >( *storage, level );
   const uint_t numLocalDoFs  = numberOfLocalDoFs< FunctionTag_T >( *storage, level );
   const uint_t numInnerDoFs = numberOfGlobalInnerDoFs< FunctionTag_T >(*storage, level);

   //local and global DoFs are the same in non parallel case:
   WALBERLA_CHECK_EQUAL( numGlobalDoFs, numLocalDoFs )
   WALBERLA_CHECK_EQUAL( expectedGlobalDoFs, numGlobalDoFs )
   WALBERLA_CHECK_EQUAL( expectedInnerDoFs, numInnerDoFs )
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( "../../data/meshes/tri_1el.msh", 2, 15, 3 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( "../../data/meshes/tri_1el.msh", 3, 45, 21 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( "../../data/meshes/tri_2el.msh", 2, 25, 9 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( "../../data/meshes/tri_2el.msh", 3, 81, 49 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( "../../data/meshes/3D/tet_1el.msh", 2, 35, 1 );
   hyteg::testFunctionProperties< hyteg::P1FunctionTag >( "../../data/meshes/3D/tet_1el.msh", 3, 165, 35 );

   return EXIT_SUCCESS;
}