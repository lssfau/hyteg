/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/trilinos/TrilinosDirectSolver.hpp"

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/ExactStencilWeights.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void trilinosSolveScalarTest( const uint_t&   solverType,
                              const uint_t&   blockPreconditionerType,
                              const uint_t&   level,
                              const MeshInfo& meshInfo,
                              const real_t&   resEps,
                              const real_t&   errEps )
{
   WALBERLA_UNUSED( blockPreconditionerType );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P2LaplaceTrilinosSolve_Domain" );

   P2Function< real_t >   x( "x", storage, level, level );
   P2Function< real_t >   x_exact( "x_exact", storage, level, level );
   P2Function< real_t >   b( "b", storage, level, level );
   P2Function< real_t >   btmp( "btmp", storage, level, level );
   P2Function< real_t >   err( "err", storage, level, level );
   P2Function< real_t >   residuum( "res", storage, level, level );
   P2Function< real_t >   nullspace( "nullspace", storage, level, level );
   P2Function< idx_t >    numerator( "numerator", storage, level, level );

   numerator.enumerate( level );

   P2ConstantLaplaceOperator A( storage, level, level );
   P2ConstantMassOperator    M( storage, level, level );

   std::function< real_t( const Point3D& ) > exact = []( const Point3D& xx ) {
      return ( 1.0 / 2.0 ) * std::sin( 2 * xx[0] ) * std::sinh( xx[1] ) * std::cos( xx[2] );
   };

   std::function< real_t( const Point3D& ) > rhs = []( const Point3D& xx ) {
      return 4 * std::sin( xx[0] ) * std::cos( xx[0] ) * std::sinh( xx[1] ) * std::cos( xx[2] );
   };

   btmp.interpolate( rhs, level, Inner );
   M.apply( btmp, b, level, All );
   x.interpolate( exact, level, DirichletBoundary );
   x_exact.interpolate( exact, level );

   //   VTKOutput vtkOutput( "../../output", "P2LaplaceTrilinosSolve", storage );
   //   vtkOutput.add( x );
   //   vtkOutput.add( x_exact );
   //   vtkOutput.add( err );
   //   vtkOutput.add( b );
   //   vtkOutput.write( level, 0 );

   uint_t localDoFs  = numberOfLocalDoFs< P2FunctionTag >( *storage, level );
   uint_t globalDoFs = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );

   WALBERLA_LOG_INFO( "localDoFs: " << localDoFs << " globalDoFs: " << globalDoFs );

   trilinos::TrilinosDirectSolver< P2ConstantLaplaceOperator > solver_0(
       trilinos::TrilinosDirectSolverType::MUMPS, storage, level, Inner | NeumannBoundary );

   walberla::WcTimer timer;
   switch ( solverType )
   {
   case 0:
      WALBERLA_LOG_INFO_ON_ROOT( "MUMPS ..." )
      solver_0.solve( A, x, b, level );
      break;
   default:
      WALBERLA_ABORT( "No solver selected" );
      break;
   }

   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   A.apply( x, btmp, level, Inner );
   residuum.assign( {1.0, -1.0}, {b, btmp}, level, Inner );

   err.assign( {1.0, -1.0}, {x, x_exact}, level );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs );
   real_t residuum_l2  = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum = " << residuum_l2 );

   //   vtkOutput.write( level, 1 );

   WALBERLA_CHECK_LESS( residuum_l2, resEps );
   WALBERLA_CHECK_LESS( discr_l2_err, errEps );

   auto tt = storage->getTimingTree()->getReduced().getCopyWithRemainder();
   // WALBERLA_LOG_INFO_ON_ROOT( tt );
}

void trilinosSolveStokesTest( const uint_t&   solverType,
                              const uint_t&   blockPreconditionerType,
                              const uint_t&   level,
                              const MeshInfo& meshInfo,
                              const real_t&   resEps,
                              const real_t&   errEpsUSum,
                              const real_t&   errEpsP )
{
   WALBERLA_UNUSED( blockPreconditionerType );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P2P1Stokes3DTrilinosSolve_Domain" );

   P2P1TaylorHoodFunction< real_t >   x( "x", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   x_exact( "x_exact", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   b( "b", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   btmp( "btmp", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   err( "err", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   residuum( "res", storage, level, level );
   P2P1TaylorHoodFunction< real_t >   nullspace( "nullspace", storage, level, level );
   P2P1TaylorHoodFunction< idx_t >    numerator( "numerator", storage, level, level );

   numerator.enumerate( level );

   P2P1TaylorHoodStokesOperator A( storage, level, level );
   P2ConstantMassOperator       M( storage, level, level );
#if 1
   std::function< real_t( const Point3D& ) > exactU = []( const Point3D& xx ) {
      return -real_c( 4 ) * std::cos( real_c( 4 ) * xx[2] );
   };
   std::function< real_t( const Point3D& ) > exactV = []( const Point3D& xx ) {
      return real_c( 8 ) * std::cos( real_c( 8 ) * xx[0] );
   };
   std::function< real_t( const Point3D& ) > exactW = []( const Point3D& xx ) {
      return -real_c( 2 ) * std::cos( real_c( 2 ) * xx[1] );
   };

   std::function< real_t( const Point3D& ) > exactP = []( const Point3D& xx ) {
      return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] );
   };

   std::function< real_t( const Point3D& ) > forceU = []( const Point3D& xx ) {
      return 4 * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ) * std::cos( 4 * xx[0] ) - 64 * std::cos( 4 * xx[2] );
   };
   std::function< real_t( const Point3D& ) > forceV = []( const Point3D& xx ) {
      return 8 * std::sin( 4 * xx[0] ) * std::sin( 2 * xx[2] ) * std::cos( 8 * xx[1] ) + 512 * std::cos( 8 * xx[0] );
   };
   std::function< real_t( const Point3D& ) > forceW = []( const Point3D& xx ) {
      return 2 * std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::cos( 2 * xx[2] ) - 8 * std::cos( 2 * xx[1] );
   };
#endif

   //   std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& xx ) { return real_c(20) * xx[0] * xx[1] * xx[1] * xx[1]; };
   //   std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& xx ) { return real_c(5) * xx[0] * xx[0] * xx[0] * xx[0] - real_c(5) * xx[1] * xx[1] * xx[1] * xx[1]; };
   //   std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& xx ) { return real_c(60) * std::pow( xx[0], 2.0 ) * xx[1] - real_c(20) * std::pow( xx[1], 3.0 ); };
   //   std::function< real_t( const hyteg::Point3D& ) > zero =   []( const hyteg::Point3D&    ) { return real_c(0); };

   btmp.uvw().interpolate( {forceU, forceV, forceW}, level, Inner );

   M.apply( btmp.uvw()[0], b.uvw()[0], level, All );
   M.apply( btmp.uvw()[1], b.uvw()[1], level, All );
   M.apply( btmp.uvw()[2], b.uvw()[2], level, All );

   x.uvw().interpolate( {exactU, exactV, exactW}, level, DirichletBoundary );

   x_exact.uvw().interpolate( {exactU, exactV, exactW}, level );
   x_exact.p().interpolate( exactP, level );

   vertexdof::projectMean( x_exact.p(), level );

   //   VTKOutput vtkOutput( "../../output", "P2P1Stokes3DTrilinosSolve", storage );
   //   vtkOutput.add( x.u );
   //   vtkOutput.add( x.v );
   //   vtkOutput.add( x.w );
   //   vtkOutput.add( x.p() );
   //   vtkOutput.add( x_exact.u );
   //   vtkOutput.add( x_exact.v );
   //   vtkOutput.add( x_exact.w );
   //   vtkOutput.add( x_exact.p() );
   //   vtkOutput.add( err.u );
   //   vtkOutput.add( err.v );
   //   vtkOutput.add( err.w );
   //   vtkOutput.add( err.p() );
   //   vtkOutput.add( b.u );
   //   vtkOutput.add( b.v );
   //   vtkOutput.add( b.w );
   //   vtkOutput.add( b.p() );
   //   vtkOutput.write( level, 0 );

   uint_t localDoFs1         = numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   uint_t globalDoFs1        = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   uint_t globalDoFsvelocity = 3 * numberOfGlobalDoFs< P2FunctionTag >( *storage, level );

   WALBERLA_LOG_INFO( "localDoFs: " << localDoFs1 << " globalDoFs: " << globalDoFs1
                                    << ", global velocity dofs: " << globalDoFsvelocity );

   trilinos::TrilinosDirectSolver< P2P1TaylorHoodStokesOperator > solver_0(
       trilinos::TrilinosDirectSolverType::MUMPS, storage, level, Inner | NeumannBoundary );

   walberla::WcTimer timer;
   switch ( solverType )
   {
   case 0:
      WALBERLA_LOG_INFO_ON_ROOT( "MUMPS ..." )
      solver_0.solve( A, x, b, level );
      break;
   default:
      WALBERLA_ABORT( "No solver selected" );
      break;
   }

   timer.end();

   vertexdof::projectMean( x.p(), level );

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   A.apply( x, residuum, level, Inner );

   err.assign( {1.0, -1.0}, {x, x_exact}, level );

   real_t discr_l2_err     = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_w = std::sqrt( err.uvw()[2].dotGlobal( err.uvw()[2], level ) / (real_t) globalDoFs1 );
   real_t discr_l2_err_1_p = std::sqrt( err.p().dotGlobal( err.p(), level ) / (real_t) globalDoFs1 );
   real_t residuum_l2_1    = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs1 );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error w = " << discr_l2_err_1_w );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1  = " << residuum_l2_1 );

   //   vtkOutput.write( level, 1 );

   WALBERLA_CHECK_LESS( residuum_l2_1, resEps );
   WALBERLA_CHECK_LESS( discr_l2_err_1_u + discr_l2_err_1_v, errEpsUSum );
   WALBERLA_CHECK_LESS( discr_l2_err_1_p, errEpsP );

   auto tt = storage->getTimingTree()->getReduced().getCopyWithRemainder();
   // WALBERLA_LOG_INFO_ON_ROOT( tt );
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   trilinosSolveScalarTest( 0, 0, 3, MeshInfo::meshCuboid( Point3D( {0, 0, 0} ), Point3D( {1, 1, 1} ), 1, 1, 1 ), 1e-14, 1e-04 );

   //   trilinosSolveStokesTest(
   //       0, 0, 2, MeshInfo::meshCuboid( Point3D( {0, 0, 0} ), Point3D( {1, 1, 1} ), 1, 1, 1 ), 2.9e-12, 0.021, 0.33 );
   //   trilinosSolveStokesTest(
   //       1, 0, 2, MeshInfo::meshCuboid( Point3D( {0, 0, 0} ), Point3D( {1, 1, 1} ), 1, 1, 1 ), 2.9e-12, 0.021, 0.33 );
   //   trilinosSolveStokesTest(
   //       2, 0, 2, MeshInfo::meshCuboid( Point3D( {0, 0, 0} ), Point3D( {1, 1, 1} ), 1, 1, 1 ), 2.9e-12, 0.021, 0.33 );
   //   trilinosSolveStokesTest(
   //       2, 1, 2, MeshInfo::meshCuboid( Point3D( {0, 0, 0} ), Point3D( {1, 1, 1} ), 1, 1, 1 ), 2.9e-12, 0.021, 0.33 );

   return EXIT_SUCCESS;
}
