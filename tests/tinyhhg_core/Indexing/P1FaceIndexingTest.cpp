
#include <tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp>
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"

#include "core/mpi/all.h"

namespace hhg {

using indexing::Index;

void checkIndices(uint_t col, uint_t row, std::vector<uint_t> ref, uint_t type){
  std::vector<size_t> result;
  switch(type){
    //vertex
    case 0:
      for( const auto & n : hhg::vertexdof::macroface::neighborsWithCenter )
      {
        result.push_back( hhg::vertexdof::macroface::indexFromVertex<3>(col, row, n) );
      }
      break;
    case 1:
      for( const auto & n : hhg::vertexdof::macroface::neighborsFromGrayFace )
      {
        result.push_back( hhg::vertexdof::macroface::indexFromGrayFace<3>(col, row, n) );
      }
      break;
    case 2:
      for( const auto & n : hhg::vertexdof::macroface::neighborsFromBlueFace )
      {
        result.push_back( hhg::vertexdof::macroface::indexFromBlueFace<3>(col, row, n) );
      }
      break;
    default:
      WALBERLA_ABORT("wrong type")
  }

  for(size_t i = 0; i < ref.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(ref[i],result[i],"col: " << col << " row: " << row <<" i: " << i);
  }
}

void testNeighborhood2( const Index & index, const std::map< stencilDirection, uint_t > & dirToLinearIndex )
{
  const uint_t level = 2;
  for ( const auto & it : dirToLinearIndex )
  {
    const auto dir         = it.first;
    const auto linearIndex = it.second;
    WALBERLA_CHECK_EQUAL( vertexdof::macroface::indexFromVertex< level >( index.x(), index.y(), dir ), linearIndex );
  }
}

void testNeighborhood3( const Index & index, const std::map< stencilDirection, uint_t > & dirToLinearIndex )
{
  const uint_t level = 3;
  for ( const auto & it : dirToLinearIndex )
  {
    const auto dir         = it.first;
    const auto linearIndex = it.second;
    WALBERLA_CHECK_EQUAL( vertexdof::macroface::indexFromVertex< level >( index.x(), index.y(), dir ), linearIndex );
  }
}

}

//const Dir neighbors_with_center[] = {S, SE, W, C, E, NW, N};
int main(int argc, char* argv[])
{
  typedef hhg::stencilDirection sd;

  //this test is written for level 3
  walberla::mpi::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  /// CHECK VERTEX INDEXING ///
  std::vector<size_t> refOneOne = {10,1,2,11,18,17,9};
  hhg::checkIndices(1,1,refOneOne,0);
  std::vector<size_t> refFiveTwo = {22,14,15,23,29,28,21};
  hhg::checkIndices(5,2,refFiveTwo,0);

  /// CHECK CELL GRAY INDEXING ///
  hhg::checkIndices(2,2,{19,20,26},1);
  hhg::checkIndices(2,4,{32,33,37},1);

  /// CHECK CELL BLUE INDEXING ///
  hhg::checkIndices(3,3,{28,33,34},2);
  hhg::checkIndices(4,1,{14,21,22},2);

  std::map< sd, uint_t > level2FirstInner;

  level2FirstInner[ sd::VERTEX_C ]   =  6;
  level2FirstInner[ sd::VERTEX_E ]   =  7;
  level2FirstInner[ sd::VERTEX_W ]   =  5;
  level2FirstInner[ sd::VERTEX_N ]   = 10;
  level2FirstInner[ sd::VERTEX_S ]   =  1;

  level2FirstInner[ sd::VERTEX_NW ]  =  9;
  level2FirstInner[ sd::VERTEX_SE ]  =  2;

  level2FirstInner[ sd::VERTEX_BC ]  = 20;
  level2FirstInner[ sd::VERTEX_BW ]  = 19;
  level2FirstInner[ sd::VERTEX_BS ]  = 16;
  level2FirstInner[ sd::VERTEX_BSW ] = 15;

  level2FirstInner[ sd::VERTEX_FC ]  = 25;
  level2FirstInner[ sd::VERTEX_FE ]  = 26;
  level2FirstInner[ sd::VERTEX_FN ]  = 29;
  level2FirstInner[ sd::VERTEX_FNE ] = 30;

  hhg::testNeighborhood2( hhg::indexing::Index( 1, 1, 1 ), level2FirstInner );

  std::map< sd, uint_t > level2SecondInner;

  level2SecondInner[ sd::VERTEX_C ]   =  7;
  level2SecondInner[ sd::VERTEX_E ]   =  8;
  level2SecondInner[ sd::VERTEX_W ]   =  6;
  level2SecondInner[ sd::VERTEX_N ]   = 11;
  level2SecondInner[ sd::VERTEX_S ]   =  2;

  level2SecondInner[ sd::VERTEX_NW ]  = 10;
  level2SecondInner[ sd::VERTEX_SE ]  =  3;

  level2SecondInner[ sd::VERTEX_BC ]  = 21;
  level2SecondInner[ sd::VERTEX_BW ]  = 20;
  level2SecondInner[ sd::VERTEX_BS ]  = 17;
  level2SecondInner[ sd::VERTEX_BSW ] = 16;

  level2SecondInner[ sd::VERTEX_FC ]  = 26;
  level2SecondInner[ sd::VERTEX_FE ]  = 27;
  level2SecondInner[ sd::VERTEX_FN ]  = 30;
  level2SecondInner[ sd::VERTEX_FNE ] = 31;

  hhg::testNeighborhood2( hhg::indexing::Index( 2, 1, 1 ), level2SecondInner );

  std::map< sd, uint_t > level3FirstInner;

  level3FirstInner[ sd::VERTEX_C ]   = 10;
  level3FirstInner[ sd::VERTEX_E ]   = 11;
  level3FirstInner[ sd::VERTEX_W ]   =  9;
  level3FirstInner[ sd::VERTEX_N ]   = 18;
  level3FirstInner[ sd::VERTEX_S ]   =  1;

  level3FirstInner[ sd::VERTEX_NW ]  = 17;
  level3FirstInner[ sd::VERTEX_SE ]  =  2;

  level3FirstInner[ sd::VERTEX_BC ]  = 54;
  level3FirstInner[ sd::VERTEX_BW ]  = 53;
  level3FirstInner[ sd::VERTEX_BS ]  = 46;
  level3FirstInner[ sd::VERTEX_BSW ] = 45;

  level3FirstInner[ sd::VERTEX_FC ]  = 81;
  level3FirstInner[ sd::VERTEX_FE ]  = 82;
  level3FirstInner[ sd::VERTEX_FN ]  = 89;
  level3FirstInner[ sd::VERTEX_FNE ] = 90;

  hhg::testNeighborhood3( hhg::indexing::Index( 1, 1, 1 ), level3FirstInner );

}
