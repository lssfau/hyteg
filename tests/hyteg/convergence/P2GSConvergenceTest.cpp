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
#include "core/timing/Timer.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Random.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "core/Format.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hyteg;

/// smoother == 0: SOR, w == 1.0
/// smoother == 1: SOR_BW, w == 1.0
static void test( const std::string & meshFile, const uint_t & level, const uint_t & maxiter, const uint_t & smoother )
{
  WALBERLA_LOG_INFO_ON_ROOT( "Mesh: " << meshFile << ", level: " << level << ", smoother: " << (smoother == 0 ? "forward SOR" : "backward SOR") )
  const bool writeVTK = false;
  const bool printTiming = true;

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

  const auto meshInfo = MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  auto storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

  hyteg::P2ConstantLaplaceOperator L(storage, level, level);

  hyteg::P2Function< real_t > residuum  ("residuum",  storage, level, level);
  hyteg::P2Function< real_t > rhs       ("rhs",       storage, level, level);
  hyteg::P2Function< real_t > p2function("p1Function",storage, level, level);
  hyteg::P2Function< real_t > Lu        ("Lu",        storage, level, level);
  hyteg::P2Function< real_t > p2Exact   ("p1Exact",   storage, level, level);
  hyteg::P2Function< real_t > error     ("error",     storage, level, level);
  hyteg::P2Function< real_t > helperFun ("helperFun", storage, level, level);

  VTKOutput vtkOutput("../../output", "gs_P2", storage);
  vtkOutput.add( p2function );
  vtkOutput.add( p2Exact );
  vtkOutput.add( rhs );
  vtkOutput.add( residuum );
  vtkOutput.add( error );
  vtkOutput.add( helperFun );

  std::function<real_t(const hyteg::Point3D&)> exactFunction = [](const hyteg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hyteg::Point3D&)> ones  = [](const hyteg::Point3D&) { return 1.0; };
  walberla::math::seedRandomGenerator(0);
  std::function<real_t(const Point3D &)> rand = [](const Point3D &) { return walberla::math::realRandom(0.0, 1.0); };

  p2function.interpolate(exactFunction, level, hyteg::DirichletBoundary);
  p2function.interpolate(rand, level, hyteg::Inner);
  p2Exact.interpolate(exactFunction, level);

  real_t begin_res   = real_c( 0 );
  real_t abs_res_old = real_c( 0 );
  real_t rel_res     = real_c( 0 );
  real_t abs_res     = real_c( 0 );

  WALBERLA_LOG_INFO_ON_ROOT( walberla::format("%6s|%10s|%10s|%10s","iter","abs_res","rel_res","conv"));

  L.apply(p2function, Lu, level, hyteg::Inner);
  residuum.assign({1.0, -1.0}, { rhs, Lu }, level, hyteg::Inner);
  begin_res = std::sqrt(residuum.dotGlobal(residuum, level, hyteg::Inner));
  abs_res_old = begin_res;

  WALBERLA_LOG_INFO_ON_ROOT( walberla::format("%6d|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res/abs_res_old))
  walberla::WcTimer timer;

  for(uint_t i = 0; i < maxiter; ++i)
  {
    if ( writeVTK )
    {
      vtkOutput.write( level, i );
    }
    if ( smoother == 0 )
      L.smooth_sor(p2function, rhs, 1.0, level, hyteg::Inner);
    else if ( smoother == 1 )
      L.smooth_sor_backwards(p2function, rhs, 1.0, level, hyteg::Inner);

    L.apply(p2function, Lu, level, hyteg::Inner);
    residuum.assign({1.0, -1.0}, { rhs, Lu }, level, hyteg::Inner);
    abs_res = std::sqrt(residuum.dotGlobal(residuum, level, hyteg::Inner));
    rel_res = abs_res / begin_res;
    WALBERLA_LOG_INFO_ON_ROOT( walberla::format("%6d|%10.3e|%10.3e|%10.3e", i+1, abs_res, rel_res, abs_res/abs_res_old))
    WALBERLA_CHECK_LESS(abs_res,abs_res_old);
    abs_res_old = abs_res;
  }
  timer.end();

  WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
  error.assign({1.0, -1.0}, {p2function, p2Exact}, level);

  helperFun.interpolate(ones, level);
  real_t npoints = helperFun.dotGlobal(helperFun, level);

  real_t discr_l2_err = std::sqrt(error.dotGlobal(error, level) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

  if ( printTiming ) {
    walberla::WcTimingTree tt = timingTree->getReduced();
    WALBERLA_LOG_INFO_ON_ROOT(tt);
  }
}

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  test( "../../data/meshes/quad_8el.msh", 4, 20, 0 );

  test( "../../data/meshes/3D/tet_1el.msh", 4, 20, 0 );
  test( "../../data/meshes/3D/tet_1el.msh", 4, 20, 1 );
  test( "../../data/meshes/3D/regular_octahedron_8el.msh", 4, 20, 0 );
  test( "../../data/meshes/3D/regular_octahedron_8el.msh", 4, 20, 1 );


  return 0;
}
