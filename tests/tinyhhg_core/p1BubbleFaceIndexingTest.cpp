#include "tinyhhg_core/tinyhhg.hpp"

//using hhg::P1BubbleFace::index;
using namespace hhg::P1BubbleFace;

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
  std::vector<size_t> refOneOne = {10,1,2,11,18,17,9,46,54,53,82,88,81};
  std::vector<size_t> refFiveTwo = {37,32,33,38,41,40,36,73,77,76,105,107,104};
  std::vector<size_t> result;
  for(auto n : CoordsVertex::neighbors_with_center)
  {
    size_t idx = CoordsVertex::index<3>(1, 1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refOneOne.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refOneOne[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : CoordsVertex::neighbors_with_center)
  {
    size_t idx = CoordsVertex::index<3>(5, 2, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refFiveTwo.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refFiveTwo[i],result[i],"i: " << i);
  }
  result.clear();
  /// CHECK CELL GRAY ///
  std::vector<size_t> refGrayOneTwo = {61,19,25,18};
  std::vector<size_t> refGraySevenZero = {52,8,16,7};
  for(auto n : CoordsCellGray::neighbors_with_center)
  {
    size_t idx = CoordsCellGray::index<3>(2, 1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(idx);
  }
  for(size_t i = 0; i < refGrayOneTwo.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refGrayOneTwo[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : CoordsCellGray::neighbors_with_center)
  {
    size_t idx = CoordsCellGray::index<3>(0, 7, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(idx);
  }
  for(size_t i = 0; i < refGraySevenZero.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refGraySevenZero[i],result[i],"i: " << i);
  }
  result.clear();
  /// CHECK CELL BLUE ///
  std::vector<size_t> refBlueOneFour = {104,32,37,36};
  std::vector<size_t> refBlueThreeZero = {84,4,13,12};
  for(auto n : CoordsCellBlue::neighbors_with_center)
  {
    size_t idx = CoordsCellBlue::index<3>(4, 1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(idx);
  }
  for(size_t i = 0; i < refBlueOneFour.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refBlueOneFour[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : CoordsCellBlue::neighbors_with_center)
  {
    size_t idx = CoordsCellBlue::index<3>(0, 3, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(idx);
  }
  for(size_t i = 0; i < refBlueThreeZero.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refBlueThreeZero[i],result[i],"i: " << i);
  }


}
