/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Nils Kohl.
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
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/composites/P1P0StokesFunction.hpp"
#include "hyteg/composites/P1P0StokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void petscSolveTest( const uint_t&   level,
                     const MeshInfo& meshInfo,
                     const real_t&   resEps,
                     const real_t&   errEpsUSum,
                     const real_t&   errEpsP )
{
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );
   writeDomainPartitioningVTK( storage, "../../output", "P1P0Stokes2DPetscSolve_Domain" );

   P1P0StokesFunction< real_t > x( "x", storage, level, level );
   P1P0StokesFunction< real_t > x_exact( "x_exact", storage, level, level );
   P1P0StokesFunction< real_t > b( "b", storage, level, level );
   P1P0StokesFunction< real_t > err( "err", storage, level, level );
   P1P0StokesFunction< real_t > residuum( "res", storage, level, level );

   P1P0StokesOperator A( storage, level, level, 0.1 );

   // output matrix
   std::string                             fileName = "../../output/p1p0stokes.m";
   PETScSparseMatrix< P1P0StokesOperator > mat;
   P1P0StokesFunction< idx_t >             numeratorSrc( "numerator", storage, level, level );
   P1P0StokesFunction< idx_t >             numeratorDst( "numerator", storage, level, level );
   numeratorSrc.enumerate( level );
   numeratorDst.enumerate( level );
   mat.createMatrixFromOperator( A, level, numeratorSrc, numeratorDst );
   mat.print( fileName, false, PETSC_VIEWER_ASCII_MATLAB );

   std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& xx ) {
      return real_c( 20 ) * xx[0] * std::pow( xx[1], 3.0 );
   };
   std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& xx ) {
      return real_c( 5 ) * std::pow( xx[0], 4.0 ) - real_c( 5 ) * std::pow( xx[1], 4.0 );
   };
   std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& xx ) {
      return real_c( 60 ) * std::pow( xx[0], 2.0 ) * xx[1] - real_c( 20 ) * std::pow( xx[1], 3.0 );
   };
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return real_c( 0 ); };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return real_c( 1 ); };

   b.uvw().interpolate( { exactU, exactV }, level, hyteg::DirichletBoundary );
   x.uvw().interpolate( { exactU, exactV }, level, DirichletBoundary );
   x_exact.uvw().interpolate( { exactU, exactV }, level );
   x_exact.p().interpolate( exactP, level );

   VTKOutput vtkOutput( "../../output", "P1P0Stokes2DPetscSolve", storage );
   vtkOutput.add( x.uvw() );
   vtkOutput.add( x.p() );
   vtkOutput.add( x_exact.uvw() );
   vtkOutput.add( x_exact.p() );
   vtkOutput.add( err.uvw() );
   vtkOutput.add( err.p() );
   vtkOutput.add( b.uvw() );
   vtkOutput.add( b.p() );
   vtkOutput.write( level, 0 );

   //   uint_t localDoFs  = hyteg::numberOfLocalDoFs< P1P0StokesFunctionTag >( *storage, level );
   //   uint_t globalDoFs = hyteg::numberOfGlobalDoFs< P1P0StokesFunctionTag >( *storage, level );
   //
   //   WALBERLA_LOG_INFO( "local DoFs: " << localDoFs << " global DoFs: " << globalDoFs );

   PETScMinResSolver< P1P0StokesOperator > solver( storage, level, 1e-8, 1e-8, 5000 );

   walberla::WcTimer timer;
   solver.solve( A, x, b, level );
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );

   // A.apply( x, residuum, level, hyteg::Inner );

   err.assign( { 1.0, -1.0 }, { x, x_exact }, level );

   vtkOutput.write( level, 1 );

#if 0

   real_t discr_l2_err_1_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_p = std::sqrt( err.p().dotGlobal( err.p(), level ) / (real_t) globalDoFs1 );
   real_t residuum_l2_1    = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs1 );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

   WALBERLA_CHECK_LESS( residuum_l2_1, resEps );
   WALBERLA_CHECK_LESS( discr_l2_err_1_u + discr_l2_err_1_v, errEpsUSum );
   WALBERLA_CHECK_LESS( discr_l2_err_1_p, errEpsP );
#endif
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   PETScManager petscManager( &argc, &argv );

   petscSolveTest( 4,
                   hyteg::MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), hyteg::MeshInfo::CRISS, 2, 2 ),
                   1.7e-13,
                   0.025,
                   0.366 );

   return EXIT_SUCCESS;
}
