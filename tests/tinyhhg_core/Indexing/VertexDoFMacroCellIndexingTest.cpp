
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/indexing/Common.hpp"

#include "core/mpi/all.h"
#include "core/debug/CheckFunctions.h"

namespace hhg {

using indexing::Index;

void testNeighborhood2( const Index & index, const std::map< stencilDirection, uint_t > & dirToLinearIndex )
{
  const uint_t level = 2;
  for ( const auto & it : dirToLinearIndex )
  {
    const auto dir         = it.first;
    const auto linearIndex = it.second;
    WALBERLA_CHECK_EQUAL( vertexdof::macrocell::indexFromVertex( level, index.x(), index.y(), index.z(), dir ), linearIndex );
  }
}

void testNeighborhood3( const Index & index, const std::map< stencilDirection, uint_t > & dirToLinearIndex )
{
  const uint_t level = 3;
  for ( const auto & it : dirToLinearIndex )
  {
    const auto dir         = it.first;
    const auto linearIndex = it.second;
    WALBERLA_CHECK_EQUAL( vertexdof::macrocell::indexFromVertex( level, index.x(), index.y(), index.z(), dir ), linearIndex );
  }
}

static void testVertexDofMacroCellIndexing()
{
  typedef stencilDirection sd;

  std::map< sd, uint_t > level2Center;

  level2Center[ sd::VERTEX_C ]   = 20;
  level2Center[ sd::VERTEX_E ]   = 21;
  level2Center[ sd::VERTEX_W ]   = 19;
  level2Center[ sd::VERTEX_N ]   = 23;
  level2Center[ sd::VERTEX_S ]   = 16;

  level2Center[ sd::VERTEX_NW ]  = 22;
  level2Center[ sd::VERTEX_SE ]  = 17;

  level2Center[ sd::VERTEX_TC ]  = 29;
  level2Center[ sd::VERTEX_TW ]  = 28;
  level2Center[ sd::VERTEX_TS ]  = 26;
  level2Center[ sd::VERTEX_TSE ] = 25;

  level2Center[ sd::VERTEX_BC ]  = 6;
  level2Center[ sd::VERTEX_BE ]  = 7;
  level2Center[ sd::VERTEX_BN ]  = 10;
  level2Center[ sd::VERTEX_BNW ] = 11;

  testNeighborhood2( Index( 1, 1, 1 ), level2Center );

  std::map< sd, uint_t > level3FirstInner;

  level3FirstInner[ sd::VERTEX_C ]   = 54;
  level3FirstInner[ sd::VERTEX_E ]   = 55;
  level3FirstInner[ sd::VERTEX_W ]   = 53;
  level3FirstInner[ sd::VERTEX_N ]   = 61;
  level3FirstInner[ sd::VERTEX_S ]   = 46;

  level3FirstInner[ sd::VERTEX_NW ]  = 60;
  level3FirstInner[ sd::VERTEX_SE ]  = 47;

  level3FirstInner[ sd::VERTEX_TC ]  = 89;
  level3FirstInner[ sd::VERTEX_TW ]  = 88;
  level3FirstInner[ sd::VERTEX_TS ]  = 82;
  level3FirstInner[ sd::VERTEX_TSE ] = 81;

  level3FirstInner[ sd::VERTEX_BC ]  = 10;
  level3FirstInner[ sd::VERTEX_BE ]  = 11;
  level3FirstInner[ sd::VERTEX_BN ]  = 18;
  level3FirstInner[ sd::VERTEX_BNW ] = 19;

  testNeighborhood3( Index( 1, 1, 1 ), level3FirstInner );

  std::map< sd, uint_t > level3TopInner;

  level3TopInner[ sd::VERTEX_C ]   = 76;
  level3TopInner[ sd::VERTEX_E ]   = 77;
  level3TopInner[ sd::VERTEX_W ]   = 75;
  level3TopInner[ sd::VERTEX_N ]   = 79;
  level3TopInner[ sd::VERTEX_S ]   = 72;

  level3TopInner[ sd::VERTEX_NW ]  = 78;
  level3TopInner[ sd::VERTEX_SE ]  = 73;

  level3TopInner[ sd::VERTEX_TC ]  = 107;
  level3TopInner[ sd::VERTEX_TW ]  = 106;
  level3TopInner[ sd::VERTEX_TS ]  = 104;
  level3TopInner[ sd::VERTEX_TSE ] = 103;

  level3TopInner[ sd::VERTEX_BC ]  = 36;
  level3TopInner[ sd::VERTEX_BE ]  = 37;
  level3TopInner[ sd::VERTEX_BN ]  = 40;
  level3TopInner[ sd::VERTEX_BNW ] = 41;

  testNeighborhood3( Index( 1, 5, 1 ), level3TopInner );
}

}

int main(int argc, char* argv[])
{
  walberla::mpi::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  hhg::testVertexDofMacroCellIndexing();
}

