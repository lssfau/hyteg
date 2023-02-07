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
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void P2CGTest( const std::string& meshFile, const uint_t level, const real_t targetError, const bool localMPI )
{
   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   const auto globalNumberOfPrimitives = setupStorage.getNumberOfPrimitives();
   const auto rank                     = uint_c( walberla::mpi::MPIManager::instance()->rank() );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   // start allocation and assembly only on root
   loadbalancing::allPrimitivesOnRoot( setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P2CGConvergenceTest_domain_before_migration" );

   std::vector< PrimitiveID > localPrimitiveIDs;
   storage->getPrimitiveIDs( localPrimitiveIDs );
   if ( rank == 0 )
   {
      WALBERLA_CHECK( localPrimitiveIDs.size() == globalNumberOfPrimitives );
   }
   else
   {
      WALBERLA_CHECK( localPrimitiveIDs.empty() );
   }

   hyteg::P2ConstantLaplaceOperator L( storage, level, level );

   hyteg::P2Function< real_t > r( "r", storage, level, level );
   hyteg::P2Function< real_t > f( "f", storage, level, level );
   hyteg::P2Function< real_t > u( "u", storage, level, level );
   hyteg::P2Function< real_t > u_exact( "u_exact", storage, level, level );
   hyteg::P2Function< real_t > err( "err", storage, level, level );
   hyteg::P2Function< real_t > npoints_helper( "npoints_helper", storage, level, level );

   if ( localMPI )
   {
      u.setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );
   }

//   hyteg::VTKOutput vtkOutput( "../../output", "P2CGConvergenceTest", storage );
//   vtkOutput.add( u );
//   vtkOutput.add( u_exact );
//   vtkOutput.add( f );
//   vtkOutput.add( r );
//   vtkOutput.add( err );
//   vtkOutput.add( npoints_helper );
//   vtkOutput.write( level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > rhs   = []( const hyteg::Point3D& ) { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   auto solver = hyteg::CGSolver< hyteg::P2ConstantLaplaceOperator >( storage, level, level );

   // instead of solving directly, we migrate the primitives in parallel
   loadbalancing::distributed::roundRobin( *storage );
   writeDomainPartitioningVTK( storage, "../../output", "P2CGConvergenceTest_domain_after_migration" );
   storage->getPrimitiveIDs( localPrimitiveIDs );
   WALBERLA_CHECK( !localPrimitiveIDs.empty() );

//   vtkOutput.write( level );

   solver.solve( L, u, f, level );

   err.assign( {1.0, -1.0}, {u, u_exact}, level );
   npoints_helper.interpolate( ones, level );

   const real_t npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   const real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / npoints );

//   vtkOutput.write( level );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   WALBERLA_CHECK_LESS( discr_l2_err, targetError );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::P2CGTest( "../../data/meshes/tri_1el.msh", 3, real_c( 1e-7 ), false );
   hyteg::P2CGTest( "../../data/meshes/quad_4el.msh", 3, real_c( 4e-7 ), false );
   hyteg::P2CGTest( "../../data/meshes/annulus_coarse.msh", 2, real_c( 4e-6 ), false );
   hyteg::P2CGTest( "../../data/meshes/3D/tet_1el.msh", 2, real_c( 4e-6 ), true );
   hyteg::P2CGTest( "../../data/meshes/3D/pyramid_2el.msh", 2, real_c( 3e-5 ), false );
   hyteg::P2CGTest( "../../data/meshes/3D/regular_octahedron_8el.msh", 2, real_c( 1.7e-5 ), true );
}