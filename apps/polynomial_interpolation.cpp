#include <tinyhhg_core/tinyhhg.hpp>
#include <core/Environment.h>

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{
  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  const uint_t PolyDegree = 7;
  const uint_t PolyInterpolationLevel = 3;

  typedef Polynomial2D<PolyDegree, PolyInterpolationLevel> Polynomial;
  typedef LSQInterpolator<PolyDegree, PolyInterpolationLevel> Interpolator;

  Polynomial poly;
  Point2D xtest{{1.0/3.0, 1.0/3.0}};

  WALBERLA_LOG_INFO("eval = " << poly.eval(xtest));

  Interpolator interpolator;
  std::vector<real_t> values(Interpolator::NumVertices);

  Point2D x;
  uint_t rowsize = levelinfo::num_microvertices_per_edge(PolyInterpolationLevel);
  uint_t inner_rowsize = rowsize;
  real_t h = 1.0 / (rowsize-1);

  uint_t offset = 0;
  for (uint_t i = 0; i < rowsize; ++i) {
    x[0] = i * h;

    for (uint_t j = 0; j < inner_rowsize; ++j) {
      x[1] = j * h;

      values[offset] = 7.0; // + std::pow(3.0 * x[0], 4) + std::pow(3.0 * x[1], 4);
      ++offset;
    }

    --inner_rowsize;
  }

  interpolator.interpolate(values, poly);

  WALBERLA_LOG_INFO("eval = " << poly.eval(xtest));

  return 0;
}
