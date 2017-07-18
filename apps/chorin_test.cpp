#include <tinyhhg_core/tinyhhg.hpp>

#include <fmt/format.h>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  WALBERLA_LOG_INFO_ON_ROOT("Navier-Stokes solver");

  hhg::Mesh mesh("../data/meshes/bfs_12el_neumann.msh");

  real_t viscosity = 1e-1;
  real_t dt = 1e-4;

  size_t minLevel = 2;
  size_t maxLevel = 4;

  std::function<real_t(const hhg::Point3D&)> bc_x = [](const hhg::Point3D& x) {
    if (x[0] < 1e-8)
    {
      return 16.0 * (x[1]-0.5) * (1.0 - x[1]);
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

  hhg::P1FunctionOld u("u", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld v("v", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld p("p", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld p_rhs("p_rhs", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld p_res("p_res", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld tmp("tmp", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld tmp2("tmp2", mesh, minLevel, maxLevel);
  hhg::P1FunctionOld res("res", mesh, minLevel, maxLevel);

  hhg::P1LaplaceOperator A(mesh, minLevel, maxLevel);
  hhg::P1DivxOperator div_x(mesh, minLevel, maxLevel);
  hhg::P1DivyOperator div_y(mesh, minLevel, maxLevel);
  hhg::P1DivTxOperator divT_x(mesh, minLevel, maxLevel);
  hhg::P1DivTyOperator divT_y(mesh, minLevel, maxLevel);
  hhg::P1MassOperator mass(mesh, minLevel, maxLevel);

  real_t time = 0.0;
  size_t iter = 0;
  size_t max_cg_iter = 10000;

  auto mass_solver = hhg::CGSolver<hhg::P1FunctionOld, hhg::P1MassOperator>(mesh, minLevel, maxLevel);
  auto laplace_solver = hhg::CGSolver<hhg::P1FunctionOld, hhg::P1LaplaceOperator>(mesh, minLevel, maxLevel);

  u.interpolate(bc_x, maxLevel, hhg::DirichletBoundary);
  v.interpolate(bc_y, maxLevel, hhg::DirichletBoundary);
  p.interpolate(zero, maxLevel-1, hhg::NeumannBoundary);

  hhg::VTKWriter({&u, &v, &p, &p_rhs}, maxLevel, "../output", fmt::format("test_{:0>4}", iter));
  ++iter;

  while (time < 1.0)
  {
    u.interpolate(bc_x, maxLevel, hhg::DirichletBoundary);
    v.interpolate(bc_y, maxLevel, hhg::DirichletBoundary);

    fmt::print("predict u\n");
    A.apply(u, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp2.interpolate(zero, maxLevel);
    mass_solver.solve(mass, tmp2, tmp, res, maxLevel, 1e-8, max_cg_iter, hhg::Inner | hhg::NeumannBoundary, true); // project
    u.assign({1.0, -dt * viscosity}, {&u, &tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    fmt::print("predict v\n");
    A.apply(v, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp2.interpolate(zero, maxLevel);
    mass_solver.solve(mass, tmp2, tmp, res, maxLevel, 1e-8, max_cg_iter, hhg::Inner | hhg::NeumannBoundary, true); // project
    v.assign({1.0, -dt * viscosity}, {&v, &tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    fmt::print("solve p\n");
    p.interpolate(zero, maxLevel-1, hhg::NeumannBoundary);
    div_x.apply(u, p_rhs, maxLevel, hhg::Inner | hhg::DirichletBoundary);
    div_y.apply(v, tmp, maxLevel, hhg::Inner | hhg::DirichletBoundary);

    p_rhs.assign({ real_t(1)/dt, real_t(1)/dt }, { &p_rhs, &tmp }, maxLevel, hhg::Inner | hhg::DirichletBoundary);
    p_rhs.restrict(maxLevel, hhg::Inner | hhg::DirichletBoundary);

    laplace_solver.solve(A, p, p_rhs, p_res, maxLevel-1, 1e-8, max_cg_iter, hhg::Inner | hhg::DirichletBoundary, true);

    p.prolongate(maxLevel-1, hhg::Inner | hhg::DirichletBoundary);

    fmt::print("correct u\n");
    divT_x.apply(u, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp2.interpolate(zero, maxLevel);
    mass_solver.solve(mass, tmp2, tmp, res, maxLevel, 1e-8, max_cg_iter, hhg::Inner | hhg::NeumannBoundary, true); // project
    u.assign({1.0, -dt}, {&u, &tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    fmt::print("correct v\n");
    divT_y.apply(v, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary);
    tmp2.interpolate(zero, maxLevel);
    mass_solver.solve(mass, tmp2, tmp, res, maxLevel, 1e-8, max_cg_iter, hhg::Inner | hhg::NeumannBoundary, true); // project
    v.assign({1.0, -dt}, {&v, &tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary);

    hhg::VTKWriter({&u, &v, &p, &p_rhs}, maxLevel, "../output", fmt::format("test_{:0>4}", iter));
    time += dt;
    ++iter;

    if (iter >= 100)
      break;
  }

  return 0;
}
