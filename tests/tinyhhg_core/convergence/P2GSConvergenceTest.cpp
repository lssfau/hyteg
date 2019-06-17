#include "core/timing/Timer.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Random.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/Format.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

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

  hhg::P2ConstantLaplaceOperator L(storage, level, level);

  hhg::P2Function< real_t > residuum  ("residuum",  storage, level, level);
  hhg::P2Function< real_t > rhs       ("rhs",       storage, level, level);
  hhg::P2Function< real_t > p2function("p1Function",storage, level, level);
  hhg::P2Function< real_t > Lu        ("Lu",        storage, level, level);
  hhg::P2Function< real_t > p2Exact   ("p1Exact",   storage, level, level);
  hhg::P2Function< real_t > error     ("error",     storage, level, level);
  hhg::P2Function< real_t > helperFun ("helperFun", storage, level, level);

  VTKOutput vtkOutput("../../output", "gs_P2", storage);
  vtkOutput.add( p2function );
  vtkOutput.add( p2Exact );
  vtkOutput.add( rhs );
  vtkOutput.add( residuum );
  vtkOutput.add( error );
  vtkOutput.add( helperFun );

  std::function<real_t(const hhg::Point3D&)> exactFunction = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };
  walberla::math::seedRandomGenerator(0);
  std::function<real_t(const Point3D &)> rand = [](const Point3D &) { return walberla::math::realRandom(0.0, 1.0); };

  p2function.interpolate(exactFunction, level, hhg::DirichletBoundary);
  p2function.interpolate(rand, level, hhg::Inner);
  p2Exact.interpolate(exactFunction, level);

  real_t begin_res, abs_res_old, rel_res, abs_res = 0;

  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6s|%10s|%10s|%10s","iter","abs_res","rel_res","conv"));

  L.apply(p2function, Lu, level, hhg::Inner);
  residuum.assign({1.0, -1.0}, { rhs, Lu }, level, hhg::Inner);
  begin_res = std::sqrt(residuum.dotGlobal(residuum, level, hhg::Inner));
  abs_res_old = begin_res;

  WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res/abs_res_old))
  walberla::WcTimer timer;

  for(uint_t i = 0; i < maxiter; ++i)
  {
    if ( writeVTK )
    {
      vtkOutput.write( level, i );
    }
    if ( smoother == 0 )
      L.smooth_sor(p2function, rhs, 1.0, level, hhg::Inner);
    else if ( smoother == 1 )
      L.smooth_sor_backwards(p2function, rhs, 1.0, level, hhg::Inner);

    L.apply(p2function, Lu, level, hhg::Inner);
    residuum.assign({1.0, -1.0}, { rhs, Lu }, level, hhg::Inner);
    abs_res = std::sqrt(residuum.dotGlobal(residuum, level, hhg::Inner));
    rel_res = abs_res / begin_res;
    WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e", i+1, abs_res, rel_res, abs_res/abs_res_old))
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
