/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Nils Kohl.
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
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void P1ElementwiseCGTest( const std::string& meshFile,
                          const uint_t       level,
                          const real_t       targetError,
                          const bool         localMPI,
                          bool               storedElementMatrices )
{
   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P1CGConvergenceTest_domain" );

   hyteg::P1ElementwiseLaplaceOperator L( storage, level, level );

   if ( storedElementMatrices )
   {
      L.computeAndStoreLocalElementMatrices();
   }

   hyteg::P1Function< real_t > r( "r", storage, level, level );
   hyteg::P1Function< real_t > f( "f", storage, level, level );
   hyteg::P1Function< real_t > u( "u", storage, level, level );
   hyteg::P1Function< real_t > u_exact( "u_exact", storage, level, level );
   hyteg::P1Function< real_t > err( "err", storage, level, level );
   hyteg::P1Function< real_t > npoints_helper( "npoints_helper", storage, level, level );

   if ( localMPI )
   {
      u.setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );
   }

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > rhs   = []( const hyteg::Point3D& ) { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   auto solver = hyteg::CGSolver< hyteg::P1ElementwiseLaplaceOperator >( storage, level, level );
   solver.solve( L, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, u_exact }, level );
   npoints_helper.interpolate( ones, level );

   const real_t npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   const real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / npoints );

   //   hyteg::VTKOutput vtkOutput( "../../output", "P1ElementwiseCGConvergenceTest", storage );
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

   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/tri_1el.msh" ), 0, 1, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/quad_4el.msh" ), 0, 1, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 0, 1, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 0, 1, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 0, 1, true, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/cube_24el.msh" ), 0, 1, true, false );

   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/tri_1el.msh" ), 1, 1e-12, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/quad_4el.msh" ), 1, 2.1e-04, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 1, 1e-12, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 1, 1e-12, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 1, 3.2e-03, true, false );

   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/tri_1el.msh" ), 3, 3.3e-6, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/quad_4el.msh" ), 3, 1.4e-5, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 2, 1.1e-6, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 3, 5e-7, true, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 2, 1.5e-4, false, false );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 2, 9.1e-4, true, false );

   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/tri_1el.msh" ), 0, 1, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/quad_4el.msh" ), 0, 1, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 0, 1, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 0, 1, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 0, 1, true, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/cube_24el.msh" ), 0, 1, true, true );

   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/tri_1el.msh" ), 1, 1e-12, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/quad_4el.msh" ), 1, 2.1e-04, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 1, 1e-12, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 1, 1e-12, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 1, 3.2e-03, true, true );

   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/tri_1el.msh" ), 3, 3.3e-6, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "/quad_4el.msh" ), 3, 1.4e-5, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 2, 1.1e-6, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 3, 5e-7, true, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/pyramid_2el.msh" ), 2, 1.5e-4, false, true );
   hyteg::P1ElementwiseCGTest( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 2, 9.1e-4, true, true );
}
