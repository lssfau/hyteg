#include "tinyhhg_core/p1functionspace/P1EdgeIndex.hpp"
#include "core/mpi/all.h"

//using hhg::P1BubbleFace::index;
using namespace hhg::P1Edge;
using walberla::uint_t;

int main(int argc, char* argv[])
{

  walberla::mpi::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  WALBERLA_CHECK_EQUAL(EdgeCoordsVertex::neighbors.size(),6);
  WALBERLA_CHECK_EQUAL(EdgeCoordsVertex::neighbors_with_center.size(),7);
  WALBERLA_CHECK_EQUAL(EdgeCoordsVertex::neighbors_on_edge.size(),2);
  WALBERLA_CHECK_EQUAL(EdgeCoordsVertex::neighbors_south.size(),2);
  WALBERLA_CHECK_EQUAL(EdgeCoordsVertex::neighbors_north.size(),2);


//  std::string enumStrings[] = {
//      "VERTEX_C",
//      "VERTEX_S",
//      "VERTEX_SE",
//      "VERTEX_E",
//      "VERTEX_N",
//      "VERTEX_NW",
//      "VERTEX_W",
//      "CELL_GRAY_SE",
//      "CELL_GRAY_NE",
//      "CELL_GRAY_NW",
//      "CELL_BLUE_SE",
//      "CELL_BLUE_NW",
//      "CELL_BLUE_SW"
//  };

  /// CHECK VERTEX ///
  std::vector<uint_t> refOne = {1,9,10,2,18,17,0};
  std::vector<uint_t> refFive = {5,13,14,6,22,21,4};
  std::vector<uint_t> result;
  for(auto n : EdgeCoordsVertex::neighbors_with_center)
  {
    size_t idx = EdgeCoordsVertex::index<3>(1, n);
    result.push_back(idx);
  }
  for(size_t i = 0; i < refOne.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refOne[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : EdgeCoordsVertex::neighbors_with_center)
  {
    size_t idx = EdgeCoordsVertex::index<3>(5, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refFive.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refFive[i],result[i],"i: " << i);
  }
  result.clear();
}
