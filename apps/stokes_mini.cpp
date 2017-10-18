#include <tinyhhg_core/tinyhhg.hpp>

using walberla::real_t;

int main(int argc, char* argv[])
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/quad_4el_neumann.msh";

  hhg::MeshInfo meshInfo = hhg::MeshInfo::fromGmshFile(meshFileName);
  hhg::SetupPrimitiveStorage setupStorage(meshInfo, walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

  hhg::loadbalancing::roundRobin( setupStorage );

  size_t minLevel = 2;
  size_t maxLevel = 4;
  size_t maxiter = 1000;

  std::shared_ptr<hhg::PrimitiveStorage> storage = std::make_shared<hhg::PrimitiveStorage>(setupStorage);

  hhg::MiniStokesFunction<real_t> r("r", storage, minLevel, maxLevel);
  hhg::MiniStokesFunction<real_t> f("f", storage, minLevel, maxLevel);
  hhg::MiniStokesFunction<real_t> u("u", storage, minLevel, maxLevel);

//  hhg::MiniStokesFunction numerator("numerator", storage, minLevel, maxLevel);

  hhg::MiniStokesOperator L(storage, minLevel, maxLevel);

//  size_t num = 1;
//  numerator.enumerate(maxLevel, num);
//
//  fmt::print("num = {}\n", num);
//
//  std::ofstream fileL("../output/L.txt");
//  L.save(numerator, numerator, fileL, maxLevel, hhg::All);
//  return 0;

  std::function<real_t(const hhg::Point3D&)> bc_x = [](const hhg::Point3D& x) {
    return 4.0 * (1.0-x[1]) * x[1];
  };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> zero = [](const hhg::Point3D&) { return 0.0; };
  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };

  u.u.interpolate(zero, maxLevel);
  u.u.interpolate(bc_x, maxLevel, hhg::DirichletBoundary);
  u.v.interpolate(zero, maxLevel, hhg::DirichletBoundary);


  auto solver = hhg::MinResSolver<hhg::MiniStokesFunction<real_t>, hhg::MiniStokesOperator>(storage, minLevel, maxLevel);
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
//      fmt::print("{0:.8f}  ", face0mem[hhg::P1BubbleFace::FaceCoordsCellGray::index<maxLevel>(i, j, hhg::P1BubbleFace::FaceCoordsCellGray::CELL_GRAY_C)]);
//    }
//    std::cout << std::endl;
//  }
//  std::cout << "=======================================" << std::endl;
//
//  std::cout << "=======================================" << std::endl;
//  std::cout << "Face 0 Cell Blue: " << std::endl;
//  for (size_t i = 0; i < v_perEdge-2; ++i) {
//    for (size_t j = 0; j < v_perEdge-2 - i; ++j) {
//      fmt::print("{0:.8f}  ", face0mem[hhg::P1BubbleFace::FaceCoordsCellBlue::index<maxLevel>(i, j, hhg::P1BubbleFace::FaceCoordsCellBlue::CELL_BLUE_C)]);
//    }
//    std::cout << std::endl;
//  }
//  std::cout << "=======================================" << std::endl;

  //hhg::VTKWriter<hhg::P1Function< real_t >>({ &u.u.p1, &u.v.p1, &u.p }, maxLevel, "../output", "stokes_mini_test");
  return EXIT_SUCCESS;
}
