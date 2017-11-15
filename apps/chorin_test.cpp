#include <tinyhhg_core/tinyhhg.hpp>

#include <boost/core/null_deleter.hpp>
#include <fmt/format.h>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

  timingTree->start("Global");
  std::string meshFileName = "../data/meshes/flow_around_cylinder_fine.msh";

  real_t viscosity = 1e-3;

  bool neumann = true;
  uint_t minLevel = 2;
  uint_t maxLevel = 4;

  real_t time = 0.0;
  real_t endTime = 0.01;
  uint_t iter = 0;
  uint_t max_cg_iter = 1000;

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

#ifdef WALBERLA_BUILD_WITH_PARMETIS
  loadbalancing::distributed::parmetis( *storage );
#endif

  storage->enableGlobalTiming(timingTree);

  const real_t minimalEdgeLength = hhg::MeshQuality::getMinimalEdgeLength(storage, maxLevel);
  real_t dt = 0.01 * minimalEdgeLength;

  WALBERLA_LOG_INFO_ON_ROOT("dt = " << dt);

  std::function<real_t(const hhg::Point3D&)> bc_x = [](const hhg::Point3D& x) {
    real_t U_m = 1.5;

    if (x[0] < 1e-8)
    {
      return 4.0  * U_m * x[1]* (0.41 - x[1]) / (0.41 * 0.41);
    }
    else
    {
      return 0.0;
    }
  };

  std::function<real_t(const hhg::Point3D&)> bc_y = [](const hhg::Point3D&) {
    return 0.0;
  };

  std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> one = [](const hhg::Point3D&) { return 1.0; };

  hhg::P1Function<real_t> u("u", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> v("v", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> p("p", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> p_rhs("p_rhs", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> p_res("p_res", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> tmp("tmp", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> tmp2("tmp2", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> res("res", storage, minLevel, maxLevel);
  hhg::P1Function<real_t> ones("ones", storage, minLevel, maxLevel);

  auto u_dg = std::make_shared<hhg::DGFunction<real_t>>("u_dg", storage, minLevel, maxLevel);
  auto v_dg = std::make_shared<hhg::DGFunction<real_t>>("v_dg", storage, minLevel, maxLevel);

  auto u_dg_old = std::make_shared<hhg::DGFunction<real_t>>("u_dg", storage, minLevel, maxLevel);
  auto v_dg_old = std::make_shared<hhg::DGFunction<real_t>>("v_dg", storage, minLevel, maxLevel);

  hhg::P1LaplaceOperator A(storage, minLevel, maxLevel);
  hhg::P1DivxOperator div_x(storage, minLevel, maxLevel);
  hhg::P1DivyOperator div_y(storage, minLevel, maxLevel);
  hhg::P1DivTxOperator divT_x(storage, minLevel, maxLevel);
  hhg::P1DivTyOperator divT_y(storage, minLevel, maxLevel);
  hhg::P1LumpedInvMassOperator invDiagMass(storage, minLevel, maxLevel);

  std::array<std::shared_ptr<hhg::P1Function<real_t>>, 2> velocity{{std::shared_ptr<hhg::P1Function<real_t>>(&u, boost::null_deleter()), std::shared_ptr<hhg::P1Function<real_t>>(&v, boost::null_deleter())}};
  hhg::DGUpwindOperator<hhg::P1Function<real_t>> N(storage, velocity, minLevel, maxLevel);

  typedef hhg::CGSolver<hhg::P1Function<real_t>, hhg::P1LaplaceOperator> CoarseSolver;
  auto coarseLaplaceSolver = std::make_shared<CoarseSolver>(storage, minLevel, minLevel);
  typedef GMultigridSolver<hhg::P1Function<real_t>, hhg::P1LaplaceOperator, CoarseSolver> LaplaceSover;
  LaplaceSover laplaceSolver(storage, coarseLaplaceSolver, minLevel, maxLevel);

  u.interpolate(bc_x, maxLevel, hhg::DirichletBoundary);
  v.interpolate(bc_y, maxLevel, hhg::DirichletBoundary);
  p.interpolate(zero, maxLevel-1, hhg::NeumannBoundary);
  ones.interpolate(one, maxLevel, hhg::All);

  u_dg->projectP1(u, maxLevel, hhg::All);
  v_dg->projectP1(v, maxLevel, hhg::All);

  hhg::VTKWriter<hhg::P1Function<real_t>, hhg::DGFunction<real_t>>({&u, &v, &p}, { }, maxLevel, "../output", fmt::format("test_{:0>7}", iter));
  ++iter;

  u.interpolate(bc_x, maxLevel, hhg::DirichletBoundary);
  v.interpolate(bc_y, maxLevel, hhg::DirichletBoundary);

  while (time < endTime)
  {
    WALBERLA_LOG_INFO_ON_ROOT("time = " << time);
    u_dg_old->projectP1(u, maxLevel, hhg::All);
    v_dg_old->projectP1(v, maxLevel, hhg::All);

    N.apply(*u_dg_old, *u_dg, maxLevel, hhg::Inner, Replace);
    N.apply(*v_dg_old, *v_dg, maxLevel, hhg::Inner, Replace);

    // Predict u
    A.apply(u, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp2.integrateDG(*u_dg, ones, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp.assign({1.0, viscosity}, {&tmp2, &tmp}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    tmp2.interpolate(zero, maxLevel);
    invDiagMass.apply(tmp, tmp2, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    u.assign({1.0, -dt}, {&u, &tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    // Predict v
    A.apply(v, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp2.integrateDG(*v_dg, ones, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp.assign({1.0, viscosity}, {&tmp2, &tmp}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    tmp2.interpolate(zero, maxLevel);
    invDiagMass.apply(tmp, tmp2, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    v.assign({1.0, -dt}, {&v, &tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    // Solve p
    p.interpolate(zero, maxLevel-1, hhg::NeumannBoundary);
    div_x.apply(u, p_rhs, maxLevel, hhg::Inner | hhg::DirichletBoundary);
    div_y.apply(v, tmp, maxLevel, hhg::Inner | hhg::DirichletBoundary);

    p_rhs.assign({ real_t(1)/dt, real_t(1)/dt }, { &p_rhs, &tmp }, maxLevel, hhg::Inner | hhg::DirichletBoundary);
    p_rhs.restrict(maxLevel, hhg::Inner | hhg::DirichletBoundary);

    if (!neumann) {
      hhg::projectMean(p_rhs, tmp, maxLevel-1);
    }

    laplaceSolver.solve(A, p, p_rhs, p_res, maxLevel-1, 1e-3, max_cg_iter, hhg::Inner | hhg::DirichletBoundary, LaplaceSover::CycleType::VCYCLE, false);

    if (!neumann) {
      hhg::projectMean(p, tmp, maxLevel - 1);
    }

    p.prolongate(maxLevel-1, hhg::Inner | hhg::DirichletBoundary);

    // Correct u
    divT_x.apply(p, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp2.interpolate(zero, maxLevel);
    invDiagMass.apply(tmp, tmp2, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    u.assign({1.0, -dt}, {&u, &tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    // Correct v
    divT_y.apply(p, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp2.interpolate(zero, maxLevel);
    invDiagMass.apply(tmp, tmp2, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    v.assign({1.0, -dt}, {&v, &tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    if (iter % 100 == 0) {
      hhg::VTKWriter < hhg::P1Function < real_t > , hhg::DGFunction < real_t >> ({ &u, &v, &p }, {},maxLevel,
                                                                                 "../output", fmt::format("test_{:0>7}",
                                                                                                          iter));
    }
    time += dt;
    ++iter;

    u_dg_old.swap(u_dg);
    v_dg_old.swap(v_dg);
  }

  timingTree->stop("Global");
  auto reduced_tt = timingTree->getReduced();
  WALBERLA_LOG_INFO_ON_ROOT(reduced_tt);

  return 0;
}
