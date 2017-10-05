#include "core/Environment.h"

#include "tinyhhg_core/dgfunctionspace/DGFaceIndex.hpp"
#include "tinyhhg_core/dgfunctionspace/DGEdgeIndex.hpp"

//using hhg::P1BubbleFace::index;
using namespace hhg::DGFace;

int main(int argc, char* argv[]) {
  //this test is written for level 3
  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  /// CHECK VERTEX this is the same as for bubble///
  std::vector <size_t> refOneOne = {1, 9, 8, 37, 43, 36};
  std::vector <size_t> refTwoFive = {28, 32, 31, 60, 62, 59};
  std::vector <size_t> result;
  for(auto n : hhg::BubbleFace::neighbors)
  {
    size_t idx = indexDGFaceFromVertex<3>(1, 1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refOneOne.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refOneOne[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : hhg::BubbleFace::neighbors)
  {
    size_t idx = indexDGFaceFromVertex<3>(2, 5, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refTwoFive.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refTwoFive[i],result[i],"i: " << i);
  }
  /// END CHECK VERTEX ///
  result.clear();


  WALBERLA_CHECK_EQUAL(50, indexDGfaceFromGrayDGface<3>(1, 3, hhg::stencilDirection::CELL_BLUE_S));
  WALBERLA_CHECK_EQUAL(55, indexDGfaceFromGrayDGface<3>(2, 3, hhg::stencilDirection::CELL_BLUE_W));
  WALBERLA_CHECK_EQUAL(48, indexDGfaceFromGrayDGface<3>(5, 1, hhg::stencilDirection::CELL_BLUE_E));


  WALBERLA_CHECK_EQUAL(27, indexDGFaceFromBlueDGface<3>(1, 3, hhg::stencilDirection::CELL_GRAY_N));
  WALBERLA_CHECK_EQUAL(23, indexDGFaceFromBlueDGface<3>(2, 3, hhg::stencilDirection::CELL_GRAY_W));
  WALBERLA_CHECK_EQUAL(14, indexDGFaceFromBlueDGface<3>(5, 1, hhg::stencilDirection::CELL_GRAY_E));
  WALBERLA_CHECK_EQUAL(13, indexDGFaceFromBlueDGface<3>(5, 0, hhg::stencilDirection::CELL_GRAY_N));
  WALBERLA_CHECK_EQUAL(14, indexDGFaceFromBlueDGface<3>(6, 0, hhg::stencilDirection::CELL_GRAY_N));
  WALBERLA_CHECK_EQUAL(7 , indexDGFaceFromBlueDGface<3>(6, 0, hhg::stencilDirection::CELL_GRAY_E));
  WALBERLA_CHECK_EQUAL(6 , indexDGFaceFromBlueDGface<3>(6, 0, hhg::stencilDirection::CELL_GRAY_W));

  /// EDGE INDEXING ///
  /// CHECK VERTEX ///
  std::vector<size_t> refOne = {2,17,15,0,1,16};
  std::vector<size_t> refFive = {10,25,23,8,9,24};
  for(auto n : hhg::BubbleFace::neighbors)
  {
    size_t idx = hhg::DGEdge::indexDGFaceFromVertex<3>(1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refOne.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refOne[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : hhg::BubbleFace::neighbors)
  {
    size_t idx = hhg::DGEdge::indexDGFaceFromVertex<3>(5, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refFive.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refFive[i],result[i],"i: " << i);
  }
  result.clear();

}