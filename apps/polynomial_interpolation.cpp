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
  const uint_t PolyInterpolationLevel = 4;

  typedef Polynomial2D<PolyDegree, PolyInterpolationLevel> Polynomial;
  typedef LSQInterpolator<PolyDegree, PolyInterpolationLevel> Interpolator;

  Polynomial poly;
  Point2D xtest{{1.0/3.0, 1.0/3.0}};

  WALBERLA_LOG_DEVEL("NumVertices = " << Interpolator::NumVertices);

  WALBERLA_LOG_INFO("eval = " << poly.eval(xtest));

  Interpolator interpolator;
  std::vector<real_t> values(Interpolator::NumVertices);

  Point2D x;
  uint_t rowsize = levelinfo::num_microvertices_per_edge(PolyInterpolationLevel);
  real_t h = real_c(1.0) / real_c(rowsize-1);

  uint_t offset = 0;
  for (uint_t i = 0; i < rowsize-3; ++i) {
    x[1] = i * h + h;

    for (uint_t j = 0; j < rowsize-2-i; ++j) {
      x[0] = j * h + 0.5 * h;

      values[offset] = 7.0 + std::pow(3.0 * x[0], 4) + std::pow(3.0 * x[1], 4);
      ++offset;
    }
  }

  interpolator.interpolate(values, poly);

  WALBERLA_LOG_INFO("eval = " << poly.eval(xtest));

  return 0;
}
