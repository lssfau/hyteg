#include <tinyhhg_core/tinyhhg.hpp>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  hhg::Mesh mesh("../data/meshes/tri_1el_neumann.msh");

  size_t minLevel = 2;
  const size_t maxLevel = 4;
  size_t maxiter = 1000;

  hhg::MiniStokesFunction r("r", mesh, minLevel, maxLevel);
  hhg::MiniStokesFunction f("f", mesh, minLevel, maxLevel);
  hhg::MiniStokesFunction u("u", mesh, minLevel, maxLevel);

  hhg::MiniStokesOperator L(mesh, minLevel, maxLevel);

  std::function<real_t(const hhg::Point3D&)> bc_x = [](const hhg::Point3D& x) {
    return 4.0 * (1.0-x[1]) * x[1];
  };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate(ones, maxLevel);

  size_t npoints = (size_t) u.u.dot(u.u, maxLevel);
  size_t real_npoints = hhg::levelinfo::num_microvertices_per_face(maxLevel);

  WALBERLA_LOG_INFO_ON_ROOT(fmt::format("npoints = {}\nreal_npoints = {}\nright = {}\n", npoints, real_npoints, npoints == real_npoints))

  u.u.interpolate(zero, maxLevel);
  u.u.interpolate(bc_x, maxLevel, hhg::DirichletBoundary);
  u.v.interpolate(zero, maxLevel, hhg::DirichletBoundary);

  auto solver = hhg::MinResSolver<hhg::MiniStokesFunction, hhg::MiniStokesOperator>(mesh, minLevel, maxLevel);
  solver.solve(L, u, f, r, maxLevel, 1e-12, maxiter, hhg::Inner | hhg::NeumannBoundary, true);
//
//  for (auto vertex: u.v.mesh.vertices) {
//    hhg::P1BubbleVertex::printFunctionMemory(vertex,u.v.memory_id,maxLevel);
//  }
//  for (auto vertex: u.u.mesh.vertices) {
//    hhg::P1BubbleVertex::printFunctionMemory(vertex,u.u.memory_id,maxLevel);
//  }
//  for (auto vertex: u.p.mesh.vertices) {
//    hhg::P1BubbleVertex::printFunctionMemory(vertex,u.p.memory_id,maxLevel);
//  }
//
//  size_t v_perEdge = hhg::levelinfo::num_microvertices_per_edge(maxLevel);
//  auto& face0mem = hhg::P1Bubble::getFaceFunctionMemory(mesh.faces[0], u.u.memory_id)->data[maxLevel];
//  std::cout << "=======================================" << std::endl;
//  std::cout << "Face 0 Cell Gray: " << std::endl;
//  for (size_t i = 0; i < v_perEdge-1; ++i) {
//    for (size_t j = 0; j < v_perEdge-1 - i; ++j) {
//      fmt::print("{0:.8f}  ", face0mem[hhg::P1BubbleFace::CoordsCellGray::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsCellGray::CELL_GRAY_C)]);
//    }
//    std::cout << std::endl;
//  }
//  std::cout << "=======================================" << std::endl;
//
//  std::cout << "=======================================" << std::endl;
//  std::cout << "Face 0 Cell Blue: " << std::endl;
//  for (size_t i = 0; i < v_perEdge-2; ++i) {
//    for (size_t j = 0; j < v_perEdge-2 - i; ++j) {
//      fmt::print("{0:.8f}  ", face0mem[hhg::P1BubbleFace::CoordsCellBlue::index<maxLevel>(i, j, hhg::P1BubbleFace::CoordsCellBlue::CELL_BLUE_C)]);
//    }
//    std::cout << std::endl;
//  }
//  std::cout << "=======================================" << std::endl;

  hhg::VTKWriter({ &u.u, &u.v, &u.p }, maxLevel, "../output", "stokes_mini_test");
  return EXIT_SUCCESS;
}
