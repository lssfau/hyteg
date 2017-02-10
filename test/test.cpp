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
  hhg::Mesh mesh("../data/meshes/quad_4el.msh");

  size_t minLevel = 2;
  size_t maxLevel = 3;

  hhg::P1FunctionSpace V = hhg::P1FunctionSpace(mesh);

  std::function<double(const hhg::Point3D&)> expr = [](const hhg::Point3D& x) { return 1.0; };
  hhg::Function<hhg::P1FunctionSpace> u("u", V, minLevel, maxLevel);
  hhg::Function<hhg::P1FunctionSpace> Lu("Lu", V, minLevel, maxLevel);

  u.interpolate(expr, maxLevel);
  hhg::P1LaplaceOperator L(V, minLevel, maxLevel);

  hhg::VTKWriter({&u}, maxLevel, "../output", "test");
  MPI_Finalize();
  return 0;
}