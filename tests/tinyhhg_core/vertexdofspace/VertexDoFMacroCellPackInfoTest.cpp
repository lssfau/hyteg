
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/vtkwriter.hpp"

namespace hhg {

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;

static void testVertexDoFMacroCellPackInfo()
{
  const uint_t level = 3;

  auto storage = PrimitiveStorage::createFromGmshFile( "../../data/meshes/3D/tet_1el.msh" );

  auto x = std::make_shared< P1Function< real_t > >( "x", storage, level, level );

  // Writing 1.0 to all macro-faces
  for ( const auto & f : storage->getFaces() )
  {
    auto faceData = f.second->getData( x->getFaceDataID() )->getPointer( level );
    for ( const auto & it : vertexdof::macroface::Iterator( level ) )
    {
      faceData[ vertexdof::macroface::indexFromVertex< level >( it.x(), it.y(), stencilDirection::VERTEX_C ) ] = 1.0;
    }
  }

  // Macro-cell should still be zero everywhere - particularly also at the borders
  for ( const auto & f : storage->getCells() )
  {
    auto cellData = f.second->getData( x->getCellDataID() )->getPointer( level );
    for ( const auto & it : vertexdof::macrocell::BorderIterator( level, 0, 1, 2 ) )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( cellData[ vertexdof::macrocell::indexFromVertex< level >( it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ], 0.0 );
    }
    for ( const auto & it : vertexdof::macrocell::BorderIterator( level, 1, 2, 3 ) )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( cellData[ vertexdof::macrocell::indexFromVertex< level >( it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ], 0.0 );
    }
  }

  // Communicating macro-face data to the macro-cell
  x->getCommunicator( level )->setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );
  x->getCommunicator( level )->template communicate< Face, Cell >();

  // Macro-cell DoFs are 1.0 at the borders now
  for ( const auto & f : storage->getCells() )
  {
    auto cellData = f.second->getData( x->getCellDataID() )->getPointer( level );
    for ( const auto & it : vertexdof::macrocell::BorderIterator( level, 0, 1, 2 ) )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( cellData[ vertexdof::macrocell::indexFromVertex< level >( it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ], 1.0 );
    }
    for ( const auto & it : vertexdof::macrocell::BorderIterator( level, 1, 2, 3 ) )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( cellData[ vertexdof::macrocell::indexFromVertex< level >( it.x(), it.y(), it.z(), stencilDirection::VERTEX_C ) ], 1.0 );
    }
  }

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testVertexDoFMacroCellPackInfo();

   return EXIT_SUCCESS;
}
