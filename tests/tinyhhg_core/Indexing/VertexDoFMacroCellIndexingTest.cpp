
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

  level2Center[ sd::VERTEX_BC ]  = 29;
  level2Center[ sd::VERTEX_BW ]  = 28;
  level2Center[ sd::VERTEX_BS ]  = 26;
  level2Center[ sd::VERTEX_BSW ] = 25;

  level2Center[ sd::VERTEX_FC ]  = 6;
  level2Center[ sd::VERTEX_FE ]  = 7;
  level2Center[ sd::VERTEX_FN ]  = 10;
  level2Center[ sd::VERTEX_FNE ] = 11;

  testNeighborhood2( Index( 1, 1, 1 ), level2Center );

  std::map< sd, uint_t > level3FirstInner;

  level3FirstInner[ sd::VERTEX_C ]   = 54;
  level3FirstInner[ sd::VERTEX_E ]   = 55;
  level3FirstInner[ sd::VERTEX_W ]   = 53;
  level3FirstInner[ sd::VERTEX_N ]   = 61;
  level3FirstInner[ sd::VERTEX_S ]   = 46;

  level3FirstInner[ sd::VERTEX_NW ]  = 60;
  level3FirstInner[ sd::VERTEX_SE ]  = 47;

  level3FirstInner[ sd::VERTEX_BC ]  = 89;
  level3FirstInner[ sd::VERTEX_BW ]  = 88;
  level3FirstInner[ sd::VERTEX_BS ]  = 82;
  level3FirstInner[ sd::VERTEX_BSW ] = 81;

  level3FirstInner[ sd::VERTEX_FC ]  = 10;
  level3FirstInner[ sd::VERTEX_FE ]  = 11;
  level3FirstInner[ sd::VERTEX_FN ]  = 18;
  level3FirstInner[ sd::VERTEX_FNE ] = 19;

  testNeighborhood3( Index( 1, 1, 1 ), level3FirstInner );

  std::map< sd, uint_t > level3TopInner;

  level3TopInner[ sd::VERTEX_C ]   = 76;
  level3TopInner[ sd::VERTEX_E ]   = 77;
  level3TopInner[ sd::VERTEX_W ]   = 75;
  level3TopInner[ sd::VERTEX_N ]   = 79;
  level3TopInner[ sd::VERTEX_S ]   = 72;

  level3TopInner[ sd::VERTEX_NW ]  = 78;
  level3TopInner[ sd::VERTEX_SE ]  = 73;

  level3TopInner[ sd::VERTEX_BC ]  = 107;
  level3TopInner[ sd::VERTEX_BW ]  = 106;
  level3TopInner[ sd::VERTEX_BS ]  = 104;
  level3TopInner[ sd::VERTEX_BSW ] = 103;

  level3TopInner[ sd::VERTEX_FC ]  = 36;
  level3TopInner[ sd::VERTEX_FE ]  = 37;
  level3TopInner[ sd::VERTEX_FN ]  = 40;
  level3TopInner[ sd::VERTEX_FNE ] = 41;

  testNeighborhood3( Index( 1, 5, 1 ), level3TopInner );
}

}

int main(int argc, char* argv[])
{
  walberla::mpi::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();
  hhg::testVertexDofMacroCellIndexing();
}

