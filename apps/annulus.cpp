/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Daniel Drzisga, Dominik Thoennes, Marcus Mohr.
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
#include <cfenv>

#include <core/timing/Timer.h>
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"

#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p1_polar_laplacian.h"


#include "hyteg/geometry/AnnulusMap.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;
using walberla::math::pi;

using namespace hyteg;


int main(int argc, char* argv[])
{

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // setup mesh and storage
  walberla::Environment walberlaEnv( argc, argv );
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName( "theTriangle.msh" );
  // MeshInfo meshInfo = MeshInfo::meshAnnulus( 1.0, 2.0, 0.25 * pi, 0.75 * pi, MeshInfo::CRISS, 6, 2 );
  MeshInfo meshInfo = MeshInfo::meshAnnulus( 1.0, 2.0, 0.0, 2.0 * pi, MeshInfo::CROSS, 12, 2 );

  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  // set geometry map for blending
  for( auto it : setupStorage.getFaces() ) {
    Face& face = *it.second;
    setupStorage.setGeometryMap( face.getID(), std::make_shared< AnnulusMap >( face ) );
  }

  for( auto it : setupStorage.getEdges() ) {
    Edge& edge = *it.second;
    setupStorage.setGeometryMap( edge.getID(), std::make_shared< AnnulusMap >( edge, setupStorage ) );
  }

  loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  uint_t level = 2;
  P1Function< real_t > newRadius( "new radius", storage, level, level );

  std::function<real_t(const Point3D&)> newRadiusFunc =
    []( const Point3D& x ) {
    return std::sqrt( x[0]*x[0] + x[1]*x[1] );
  };

  newRadius.interpolate( newRadiusFunc, level );

  VTKOutput vtkOutput( "../output", "annulus", storage );
  vtkOutput.add( newRadius );
  vtkOutput.write( level );

  return 0;
}
