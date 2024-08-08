/*
 * Copyright (c) 2017-2019 Nils Kohl.
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

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void PrimitiveMigrationMatMulTest( const std::string& meshFile, const uint_t level, const bool localMPI )
{
   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   const auto globalNumberOfPrimitives = setupStorage.getNumberOfPrimitives();
   const auto rank                     = uint_c( walberla::mpi::MPIManager::instance()->rank() );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   // start allocation and assembly only on root
   loadbalancing::allPrimitivesOnRoot( setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "PrimitiveMigrationMatMulTest_domain_before_migration" );

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

   P2ConstantLaplaceOperator L( storage, level, level );

   P2Function< real_t > u( "u", storage, level, level );
   P2Function< real_t > f_before_migration( "f_before_migration", storage, level, level );
   P2Function< real_t > f_after_migration( "f_after_migration", storage, level, level );
   P2Function< real_t > f_error( "f_error", storage, level, level );

   if ( localMPI )
   {
      u.setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );
   }

//   hyteg::VTKOutput vtkOutput( "../../output", "PrimitiveMigrationMatMulTest", storage );
//   vtkOutput.add( u );
//   vtkOutput.add( f_before_migration );
//   vtkOutput.add( f_after_migration );
//   vtkOutput.add( f_error );
//   vtkOutput.write( level, 0 );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };

   u.interpolate( exact, level, All );

   // perform matrix-vector multiplication before migration and calculate norm
   L.apply( u, f_before_migration, level, All );
   const auto u_norm_before = u.dotGlobal( u, level, All );
   const auto f_norm_before = f_before_migration.dotGlobal( f_before_migration, level, All );

//   vtkOutput.write( level, 1 );

   // we now migrate the primitives in parallel now
   loadbalancing::distributed::roundRobin( *storage );
   writeDomainPartitioningVTK( storage, "../../output", "PrimitiveMigrationMatMulTest_domain_after_migration" );
   storage->getPrimitiveIDs( localPrimitiveIDs );
   WALBERLA_CHECK( !localPrimitiveIDs.empty() );

   L.apply( u, f_after_migration, level, All );
   const auto u_norm_after = u.dotGlobal( u, level, All );
   const auto f_norm_after = f_after_migration.dotGlobal( f_after_migration, level, All );

   f_error.assign( {1.0, -1.0}, {f_before_migration, f_after_migration}, level, All );

//   vtkOutput.write( level, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "mesh: " << meshFile << ", norm difference u: " << u_norm_after - u_norm_before );
   WALBERLA_LOG_INFO_ON_ROOT( "mesh: " << meshFile << ", norm difference f: " << f_norm_after - f_norm_before );

   WALBERLA_CHECK_FLOAT_EQUAL( u_norm_before, u_norm_after );
   WALBERLA_CHECK_FLOAT_EQUAL( f_norm_before, f_norm_after );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PrimitiveMigrationMatMulTest( hyteg::prependHyTeGMeshDir( "2D/tri_1el.msh" ), 3, false );
   hyteg::PrimitiveMigrationMatMulTest( hyteg::prependHyTeGMeshDir( "2D/quad_4el.msh" ), 3, false );
   hyteg::PrimitiveMigrationMatMulTest( hyteg::prependHyTeGMeshDir( "2D/annulus_coarse.msh" ), 2, false );
   hyteg::PrimitiveMigrationMatMulTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 3, true );
   hyteg::PrimitiveMigrationMatMulTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 2, false );
   hyteg::PrimitiveMigrationMatMulTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 2, true );
   hyteg::PrimitiveMigrationMatMulTest( hyteg::prependHyTeGMeshDir( "3D/cube_24el.msh" ), 2, true );
}
