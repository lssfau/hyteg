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
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"

#include "constantStencilOperator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv( argc, argv );
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  const uint_t      lowerLevel       = 3;
  const uint_t      higherLevel     = lowerLevel + 1;
  const std::string meshFile        = "../../data/meshes/3D/regular_octahedron_8el.msh";
  const bool        writeVTK        = false;
  const bool        enableChecks    = true;

  const auto meshInfo = MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

  WALBERLA_CHECK( storage->hasGlobalCells() );

  writeDomainPartitioningVTK( storage, "../../output", "P1_CG_3D_convergence_partitioning" );

  P1ConstantLaplaceOperator laplaceOperator3D( storage, lowerLevel, higherLevel );

  std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D & p ) -> real_t
  {
    return sin(p[0]) * sinh(p[1]) * p[2];
  };

  std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D & ) -> real_t
  {
    return 0.0;
  };

  std::function< real_t( const hyteg::Point3D& ) > one = []( const hyteg::Point3D & ) -> real_t
  {
      return 1.0;
  };

  std::function< real_t( const hyteg::Point3D& ) > rand = []( const hyteg::Point3D & ) -> real_t
  {
    return real_c( walberla::math::realRandom( 0.0, 1.0 ) );
  };

  hyteg::P1Function< real_t > res( "r", storage, lowerLevel, higherLevel );
  hyteg::P1Function< real_t > f( "f", storage, lowerLevel, higherLevel );
  hyteg::P1Function< real_t > u( "u", storage, lowerLevel, higherLevel );
  hyteg::P1Function< real_t > uExact( "u_exact", storage, lowerLevel, higherLevel );
  hyteg::P1Function< real_t > err( "err", storage, lowerLevel, higherLevel );
  hyteg::P1Function< real_t > oneFunction( "oneFunction", storage, lowerLevel, higherLevel );

  u.interpolate( rand, lowerLevel, DoFType::Inner );
  u.interpolate( exact, lowerLevel, DoFType::DirichletBoundary );
  f.interpolate( zero, lowerLevel, DoFType::All );
  res.interpolate( zero, lowerLevel, DoFType::All );
  uExact.interpolate( exact, lowerLevel, DoFType::All );
  oneFunction.interpolate( one, lowerLevel, DoFType::All );

  u.interpolate( rand, higherLevel, DoFType::Inner );
  u.interpolate( exact, higherLevel, DoFType::DirichletBoundary );
  f.interpolate( zero, higherLevel, DoFType::All );
  res.interpolate( zero, higherLevel, DoFType::All );
  uExact.interpolate( exact, higherLevel, DoFType::All );
  oneFunction.interpolate( one, higherLevel, DoFType::All );

  auto solver = hyteg::CGSolver< P1ConstantLaplaceOperator >( storage, lowerLevel, higherLevel );

  WALBERLA_CHECK_LESS( lowerLevel, higherLevel );

  const real_t numPointsLowerLevel  = oneFunction.dotGlobal( oneFunction, lowerLevel,  DoFType::Inner );
  const real_t numPointsHigherLevel = oneFunction.dotGlobal( oneFunction, higherLevel, DoFType::Inner );

  VTKOutput vtkOutput("../../output", "P1CGConvergenceTest", storage);
  vtkOutput.add( u );
  vtkOutput.add( err );

  if ( writeVTK )
  {
    vtkOutput.write( higherLevel, 0 );
    vtkOutput.write( lowerLevel, 0 );
  }

  solver.solve( laplaceOperator3D, u, f, lowerLevel );
  solver.solve( laplaceOperator3D, u, f, higherLevel );

  err.assign( {1.0, -1.0}, {u, uExact}, lowerLevel );
  err.assign( {1.0, -1.0}, {u, uExact}, higherLevel );
  laplaceOperator3D.apply( u, res, lowerLevel,  DoFType::Inner );
  laplaceOperator3D.apply( u, res, higherLevel, DoFType::Inner );

  if ( writeVTK )
  {
    vtkOutput.write( higherLevel, 1 );
    vtkOutput.write( lowerLevel, 1 );
  }

  const real_t discrL2ErrLowerLevel  = std::sqrt( err.dotGlobal( err, lowerLevel,  DoFType::Inner ) / numPointsLowerLevel );
  const real_t discrL2ErrHigherLevel = std::sqrt( err.dotGlobal( err, higherLevel, DoFType::Inner ) / numPointsHigherLevel );

  const real_t discrL2ResLowerLevel  = std::sqrt( res.dotGlobal( res, lowerLevel,  DoFType::Inner ) / numPointsLowerLevel );
  const real_t discrL2ResHigherLevel = std::sqrt( res.dotGlobal( res, higherLevel, DoFType::Inner ) / numPointsHigherLevel );

  WALBERLA_LOG_INFO( "Residual L2 on level " << lowerLevel  << ": " << std::scientific << discrL2ResLowerLevel  << " | Error L2: " << discrL2ErrLowerLevel );
  WALBERLA_LOG_INFO( "Residual L2 on level " << higherLevel << ": " << std::scientific << discrL2ResHigherLevel << " | Error L2: " << discrL2ErrHigherLevel );

  auto dp = std::is_same< real_t, double >();
  if ( enableChecks )
  {
    WALBERLA_CHECK_LESS( discrL2ResLowerLevel, dp ? 6.9e-17 : 5e-8 );
    WALBERLA_CHECK_LESS( discrL2ResHigherLevel, dp ? 3.88e-17 : 3e-8 );

    // L2 err higher level ~ 0.25 * L2 err lower level
    WALBERLA_CHECK_LESS( discrL2ErrLowerLevel, 4.4e-04 );
    WALBERLA_CHECK_LESS( discrL2ErrHigherLevel, 1.0e-04 );
  }

  return 0;
}
