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

  const uint_t PolyDegree = 4;
  const uint_t PolyInterpolationLevel = 3;
  const uint_t level = 8;

  typedef Polynomial2D<PolyDegree, PolyInterpolationLevel> Polynomial;

  Polynomial poly;
  Point2D x{{1.0, 1.0}};

  WALBERLA_LOG_INFO("eval = " << poly.eval(x));

  return 0;
}
