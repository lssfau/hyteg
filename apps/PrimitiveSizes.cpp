
#include "core/Environment.h"
#include "core/DataTypes.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/FunctionTraits.hpp"

using walberla::uint_t;

namespace hhg {

template< typename FunctionTag_T >
void printPrimitiveSizes()
{
  WALBERLA_LOG_INFO_ON_ROOT( "level |            cell |            face |            edge |          vertex |" );
  WALBERLA_LOG_INFO_ON_ROOT( "------+-----------------+-----------------+-----------------+-----------------+" );
  for ( uint_t level = 2; level < 15; level++ )
  {
    const uint_t vertexSize = numberOfInnerDoFs< FunctionTag_T, Vertex >( level );
    const uint_t edgeSize = numberOfInnerDoFs< FunctionTag_T, Edge >( level );
    const uint_t faceSize = numberOfInnerDoFs< FunctionTag_T, Face >( level );
    const uint_t cellSize = numberOfInnerDoFs< FunctionTag_T, Cell >( level );
    WALBERLA_LOG_INFO_ON_ROOT( std::setw(5) << level << " | " << std::setw(15) << cellSize << " | " << std::setw(15) << faceSize << " | " << std::setw(15) << edgeSize << " | " << std::setw(15) << vertexSize << " |" );
  }
}

}

int main( int argc, char* argv[] )
{

  walberla::Environment walberlaEnv( argc, argv );
  walberla::MPIManager::instance()->useWorldComm();

  WALBERLA_LOG_INFO_ON_ROOT( " --- Primitive Sizes (number of INNER DoFs) --- " );
  WALBERLA_LOG_INFO_ON_ROOT( "P1:" );
  WALBERLA_LOG_INFO_ON_ROOT( "" );
  hhg::printPrimitiveSizes< hhg::P1FunctionTag >();
  WALBERLA_LOG_INFO_ON_ROOT( "" );
  WALBERLA_LOG_INFO_ON_ROOT( "P2:" );
  WALBERLA_LOG_INFO_ON_ROOT( "" );
  hhg::printPrimitiveSizes< hhg::P2FunctionTag >();
}