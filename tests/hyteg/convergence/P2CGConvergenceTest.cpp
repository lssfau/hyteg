/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
#include "hyteg/solvers/CGSolver.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void P2CGTest( const std::string& meshFile, const uint_t level, const real_t targetError, const bool localMPI )
{
   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P2CGConvergenceTest_domain" );

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

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > rhs   = []( const hyteg::Point3D& ) { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   auto solver = hyteg::CGSolver< hyteg::P2ConstantLaplaceOperator >( storage, level, level );
   solver.solve( L, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, u_exact }, level );
   npoints_helper.interpolate( ones, level );

   const real_t npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   const real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / npoints );

   //   hyteg::VTKOutput vtkOutput( "../../output", "P2CGConvergenceTest", storage );
   //   vtkOutput.add( u );
   //   vtkOutput.add( u_exact );
   //   vtkOutput.add( f );
   //   vtkOutput.add( r );
   //   vtkOutput.add( err );
   //   vtkOutput.add( npoints_helper );
   //   vtkOutput.write( level );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err << " (level " << level << ", mesh: " << meshFile << ")" );
   WALBERLA_CHECK_LESS( discr_l2_err, targetError );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "/2D/tri_1el.msh" ), 0, 1, false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "/2D/quad_4el.msh" ), 0, 1, false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 0, 1, false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 0, 1, false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 0, 1, true );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/cube_24el.msh" ), 0, 1, true );

   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "/2D/tri_1el.msh" ), 1, real_c( 1.5e-05 ), false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "/2D/quad_4el.msh" ), 1, real_c( 2e-05 ), false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 1, real_c( 3e-06 ), false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 1, real_c( 2e-04 ), false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 1, real_c( 1.5e-04 ), true );

   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "/2D/tri_1el.msh" ), 3, real_c( 1e-7 ), false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "/2D/quad_4el.msh" ), 3, real_c( 4e-7 ), false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 2, real_c( 3e-6 ), false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 3, real_c( 3e-7 ), true );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 2, real_c( 3e-5 ), false );
   hyteg::P2CGTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 2, real_c( 1.7e-5 ), true );
}
