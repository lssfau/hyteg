#include "core/Environment.h"

#include "tinyhhg_core/dgfunctionspace/DgFaceIndex.hpp"

//using hhg::P1BubbleFace::index;
using namespace hhg::DgFace;

int main(int argc, char* argv[]) {
  //this test is written for level 3
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

  /// CHECK VERTEX this is the same as for bubble///
  std::vector <size_t> refOneOne = {1, 9, 8, 37, 43, 36};
  std::vector <size_t> refTwoFive = {28, 32, 31, 60, 62, 59};
  std::vector <size_t> result;
  for(auto n : hhg::BubbleFace::neighbors)
  {
    size_t idx = indexDGcellFromVertex<3>(1, 1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refOneOne.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refOneOne[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : hhg::BubbleFace::neighbors)
  {
    size_t idx = indexDGcellFromVertex<3>(2, 5, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refTwoFive.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refTwoFive[i],result[i],"i: " << i);
  }
  /// END CHECK VERTEX ///
  result.clear();


  WALBERLA_CHECK_EQUAL(50,indexDGcellFromGrayDGCell<3>(1,3,hhg::stencilDirection::CELL_BLUE_S));
  WALBERLA_CHECK_EQUAL(55,indexDGcellFromGrayDGCell<3>(2,3,hhg::stencilDirection::CELL_BLUE_W));
  WALBERLA_CHECK_EQUAL(48,indexDGcellFromGrayDGCell<3>(5,1,hhg::stencilDirection::CELL_BLUE_E));


  WALBERLA_CHECK_EQUAL(27,indexDGcellFromBlueDGCell<3>(1,3,hhg::stencilDirection::CELL_GRAY_N));
  WALBERLA_CHECK_EQUAL(23,indexDGcellFromBlueDGCell<3>(2,3,hhg::stencilDirection::CELL_GRAY_W));
  WALBERLA_CHECK_EQUAL(14,indexDGcellFromBlueDGCell<3>(5,1,hhg::stencilDirection::CELL_GRAY_E));

}