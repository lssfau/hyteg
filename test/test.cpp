#include <tinyhhg.hpp>

#include <fmt/format.h>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int rk = hhg::Comm::get().rk;
  if (rk == 0)
  {
    fmt::printf("TinyHHG default test\n");
  }
  hhg::Mesh mesh("../data/meshes/bfs_126el.msh");

  size_t minLevel = 2;
  size_t maxLevel = 4;

  std::function<double(const hhg::Point3D&)> expr = [](const hhg::Point3D& x) { return cos(6*M_PI*x[0]); };
  hhg::P2Function u("u", mesh, minLevel, maxLevel);
  // hhg::P1Function Lu("Lu", mesh, minLevel, maxLevel);

  u.interpolate(expr, maxLevel);
  // hhg::P1LaplaceOperator L(mesh, minLevel, maxLevel);
  // hhg::P1MassOperator L_gen(mesh, minLevel, maxLevel);

  hhg::VTKWriter({&u}, maxLevel, "../output", "test");
  // hhg::P2Edge::print(mesh.edges[4],u.memory_id,maxLevel);
  // hhg::P2Vertex::print(mesh.vertices[0], u.memory_id, maxLevel);
  // hhg::P2Face::print(mesh.faces[1], u.memory_id, maxLevel);


  MPI_Finalize();
  return 0;
}
