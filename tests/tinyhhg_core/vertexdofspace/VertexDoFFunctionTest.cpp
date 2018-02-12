
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

static void testVertexDoFFunction()
{
  const uint_t level = 3;

  MeshInfo mesh  = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  writeDomainPartitioningVTK( storage, "../../output", "single_tet" );

  auto x = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "x", storage, level, level );
  auto y = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "y", storage, level, level );

  x->getCommunicator( level )->setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );
  y->getCommunicator( level )->setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );

  std::function< real_t( const hhg::Point3D & ) > expr = []( const hhg::Point3D & xx ) -> real_t { return real_c( (1.0L/2.0L)*sin(2*xx[0])*sinh(xx[1]) ) * real_c( xx[2] ); };

  x->interpolate( expr, level );

  x->getCommunicator( level )->template communicate< Vertex, Edge >();
  x->getCommunicator( level )->template communicate< Edge, Face >();
  x->getCommunicator( level )->template communicate< Face, Cell >();

  y->assign( { 13.0 }, { x.get() }, level );

  y->getCommunicator( level )->template communicate< Vertex, Edge >();
  y->getCommunicator( level )->template communicate< Edge, Face >();
  y->getCommunicator( level )->template communicate< Face, Cell >();

  for ( const auto & cellIt : storage->getCells() )
  {
    const Cell & cell = *cellIt.second;
    const auto xData = cell.getData( x->getCellDataID() )->getPointer( level );
    const auto yData = cell.getData( y->getCellDataID() )->getPointer( level );
    for ( const auto & it : vertexdof::macrocell::Iterator( level ) )
    {
      const uint_t idx = vertexdof::macrocell::indexFromVertex< level >( it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
      WALBERLA_CHECK_FLOAT_EQUAL( 13.0 * xData[ idx ], yData[ idx ] );
    }
  }

  VTKOutput vtkOutput( "../../output", "vertex_dof_macro_cell_test" );
  vtkOutput.set3D();
  vtkOutput.add( x );
  vtkOutput.add( y );
  vtkOutput.write( level );

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testVertexDoFFunction();

   return EXIT_SUCCESS;
}
