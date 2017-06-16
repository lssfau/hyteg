#include "tinyhhg_core/tinyhhg.hpp"

//using hhg::P1BubbleFace::index;
using namespace hhg::P1BubbleEdge;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

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
  std::vector<size_t> refOne = {1,9,10,2,33,32,0,19,42,40,17,18,41};
  std::vector<size_t> refFive = {5,13,14,6,37,36,4,27,50,48,25,26,49};
  std::vector<size_t> result;
  for(auto n : EdgeCoordsVertex::neighbors_with_center)
  {
    size_t idx = EdgeCoordsVertex::index<3>(1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
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
