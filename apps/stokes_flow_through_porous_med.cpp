#include <tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp>
#include <tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp>
#include <tinyhhg_core/VTKWriter.hpp>
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"

#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScVector.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

using walberla::real_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/porous_fine.msh";

  hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
  hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  const uint_t level           = 2;
  const uint_t maxIterations   = 1000;
  const real_t targetResidual  = 1e-12;
  const bool   usePetsc        = true;

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

#ifdef WALBERLA_BUILD_WITH_PARMETIS
  hyteg::loadbalancing::distributed::parmetis( *storage );
#endif

  hyteg::P2P1TaylorHoodFunction< real_t > r( "r", storage, level, level );
  hyteg::P2P1TaylorHoodFunction< real_t > f( "f", storage, level, level );
  hyteg::P2P1TaylorHoodFunction< real_t > u( "u", storage, level, level );

  hyteg::P2P1TaylorHoodStokesOperator L( storage, level, level );

  std::function< real_t( const hyteg::Point3D& ) > bc_x = []( const hyteg::Point3D& x ) {
      if( x[0] < 1e-8 )
      {
        return 4.0 * x[1] * (1.0 - x[1]);
      } else
      {
        return 0.0;
      }
  };
  std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
  std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
  std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

  u.u.interpolate( bc_x, level, hyteg::DirichletBoundary );
  u.v.interpolate( zero, level, hyteg::DirichletBoundary );

  hyteg::VTKOutput vtkOutput("../output", "stokes_porous_taylor_hood", storage);

  vtkOutput.add( r.u );
  vtkOutput.add( r.v );
  vtkOutput.add( r.p );

  vtkOutput.add( f.u );
  vtkOutput.add( f.v );
  vtkOutput.add( f.p );

  vtkOutput.add( u.u );
  vtkOutput.add( u.v );
  vtkOutput.add( u.p );

  timingTree->start( "Complete app" );

  vtkOutput.write( level, 0 );

#ifdef HHG_BUILD_WITH_PETSC
  if ( usePetsc )
  {
    PETScManager petscManager;
    f.u.interpolate(bc_x, level, hyteg::DirichletBoundary);
    PETScLUSolver< hyteg::P2P1TaylorHoodStokesOperator> solver( storage, level );
    solver.solve( L, u, f, level );
  }
  else
#else
  if ( usePetsc )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "HHG was not built with PETSc - solving with default solver now..." );
  }
#endif
  {
    auto solver =
         hyteg::MinResSolver< hyteg::P2P1TaylorHoodStokesOperator >( storage, level, level, maxIterations, targetResidual );

    solver.solve( L, u, f, level );
  }

  vtkOutput.write( level, 1 );

  timingTree->stop( "Complete app" );
  WALBERLA_LOG_INFO_ON_ROOT( *timingTree );

  return EXIT_SUCCESS;
}
