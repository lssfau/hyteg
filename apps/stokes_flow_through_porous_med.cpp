/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include <hyteg/composites/P2P1TaylorHoodFunction.hpp>
#include <hyteg/composites/P2P1TaylorHoodStokesOperator.hpp>
#include <hyteg/dataexport/VTKOutput.hpp>
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1StokesOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/solvers/MinresSolver.hpp"

#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"

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

  u.uvw.u.interpolate( bc_x, level, hyteg::DirichletBoundary );
  u.uvw.v.interpolate( zero, level, hyteg::DirichletBoundary );

  hyteg::VTKOutput vtkOutput("../output", "stokes_porous_taylor_hood", storage);

  vtkOutput.add( r.uvw.u );
  vtkOutput.add( r.uvw.v );
  vtkOutput.add( r.p );

  vtkOutput.add( f.uvw.u );
  vtkOutput.add( f.uvw.v );
  vtkOutput.add( f.p );

  vtkOutput.add( u.uvw.u );
  vtkOutput.add( u.uvw.v );
  vtkOutput.add( u.p );

  timingTree->start( "Complete app" );

  vtkOutput.write( level, 0 );

#ifdef HYTEG_BUILD_WITH_PETSC
  if ( usePetsc )
  {
    PETScManager petscManager;
    f.uvw.u.interpolate(bc_x, level, hyteg::DirichletBoundary);
    PETScLUSolver< hyteg::P2P1TaylorHoodStokesOperator> solver( storage, level );
    solver.solve( L, u, f, level );
  }
  else
#else
  if ( usePetsc )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "hyteg was not built with PETSc - solving with default solver now..." );
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
