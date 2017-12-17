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

  const uint_t DEGREE = 7;
  const uint_t level = 8;

  typedef Polynomial2DBasis<DEGREE> Basis;

  Basis basis;

  WALBERLA_LOG_INFO("before");
  for(auto& poly : basis.polys_) {
    WALBERLA_LOG_INFO("poly = " << poly);
  }

  basis.orthogonalize(level);

  WALBERLA_LOG_INFO("after");
  for(auto& poly : basis.polys_) {
    WALBERLA_LOG_INFO("poly = " << poly);
  }

  real_t start, end;

  real_t sp = 0.0;

  start = walberla::timing::getWcTime();

  for (uint_t i = 0; i < Polynomial2D<DEGREE>::getNumCoefficients(); ++i) {
    for (uint_t j = 0; j < Polynomial2D<DEGREE>::getNumCoefficients(); ++j) {
      sp += PolyMath::scalarProduct2D(basis.polys_[i], basis.polys_[j], level);
    }
  }

  end = walberla::timing::getWcTime();

  WALBERLA_LOG_INFO("sp = " << sp)
  WALBERLA_LOG_INFO("notgen = " << end-start << "s")

  sp = 0.0;
  start = walberla::timing::getWcTime();

  for (uint_t i = 0; i < Polynomial2D<DEGREE>::getNumCoefficients(); ++i) {
    for (uint_t j = 0; j < Polynomial2D<DEGREE>::getNumCoefficients(); ++j) {
      sp += PolyMath::scalarProduct2DHierarchical(i, j, level);
    }
  }

  end = walberla::timing::getWcTime();

  WALBERLA_LOG_INFO("sp = " << sp)
  WALBERLA_LOG_INFO("gen = " << end-start << "s")

  return 0;
}
