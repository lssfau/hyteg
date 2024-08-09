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
#include <core/timing/Timer.h>

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
#error "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON"
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );

   std::string meshFileName = prependHyTeGMeshDir( "2D/quad_8el.msh" );

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   const size_t level = 2;

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   hyteg::P1Function< real_t > x( "x", storage, level, level );
   hyteg::P1Function< real_t > x_exact( "x_exact", storage, level, level );
   hyteg::P1Function< real_t > err( "err", storage, level, level );

   hyteg::P1ConstantLaplaceOperator A( storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return xx[0] * xx[0] - xx[1] * xx[1] + 10;
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   x.interpolate( exact, level, hyteg::DirichletBoundary );
   x_exact.interpolate( exact, level );

   uint_t globalDoFs = hyteg::numberOfGlobalDoFs< P1FunctionTag >( *storage, level );
   WALBERLA_LOG_INFO_ON_ROOT( "Num dofs = " << uint_c( globalDoFs ) )

   PETScLUSolver< hyteg::P1ConstantLaplaceOperator > solver( storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Solving System" )
   walberla::WcTimer timer;
   solver.solve( A, x, x, level );
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   err.assign( { 1.0, -1.0 }, { x, x_exact }, level );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   WALBERLA_CHECK_LESS( discr_l2_err, 1e-14 );

   //  WALBERLA_LOG_INFO_ON_ROOT("Printing Solution")
   //  hyteg::VTKWriter< P1Function >({ x, x_exact, &err }, level, "../output", "exact_solver");

   return EXIT_SUCCESS;
}
