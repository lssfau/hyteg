#include "core/Environment.h"

#include "tinyhhg_core/tinyhhg.hpp"

//using hhg::P1BubbleFace::index;
using namespace hhg::BubbleFace;

int main(int argc, char* argv[])
{
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

  /// CHECK VERTEX ///
  std::vector<size_t> refOneOne = {1,9,8,37,43,36};
  std::vector<size_t> refTwoFive = {28,32,31,60,62,59};
  std::vector<size_t> result;
  for(auto n : CoordsVertex::neighbors)
  {
    size_t idx = CoordsVertex::index<3>(1, 1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refOneOne.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refOneOne[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : CoordsVertex::neighbors)
  {
    size_t idx = CoordsVertex::index<3>(2, 5, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refTwoFive.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refTwoFive[i],result[i],"i: " << i);
  }
  result.clear();
  /// CHECK CELL GRAY ///
  WALBERLA_CHECK_EQUAL(22,CoordsCellGray::index<3>(1,3,CoordsCellGray::CELL_GRAY_C));
  WALBERLA_CHECK_EQUAL(24,CoordsCellGray::index<3>(3,3,CoordsCellGray::CELL_GRAY_C));

  /// CHECK CELL BLUE //
  WALBERLA_CHECK_EQUAL(36,CoordsCellBlue::index<3>(0,0,CoordsCellBlue::CELL_BLUE_C));
  WALBERLA_CHECK_EQUAL(41,CoordsCellBlue::index<3>(5,0,CoordsCellBlue::CELL_BLUE_C));


//  /// CHECK VERTEX ITERATOR ///
//  std::vector<size_t> vertexface0 = {0,1,2,3,4,5,6,7,8};
//  std::vector<size_t> vertexface1 = {8,16,23,29,34,38,41,43,44};
//  std::vector<size_t> vertexface2 = {44,42,39,35,30,24,17,9,0};
//  uint_t counter = 0;
//  for(auto it = indexIterator(0,1, VERTEX, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexface0[counter++],*it,"counter: " << counter);
//  }
//  WALBERLA_CHECK_EQUAL_2(counter, 9);
//  counter = 0;
//  for(auto it = indexIterator(1,1, VERTEX, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexface1[counter++],*it,"counter: " << counter);
//  }
//  WALBERLA_CHECK_EQUAL_2(counter, 9);
//  counter = 0;
//  for(auto it = indexIterator(2,1, VERTEX, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexface2[counter++],*it,"counter: " << counter);
//  }
//  WALBERLA_CHECK_EQUAL_2(counter, 9);
//  /// CHECK REVERSE ///
//  counter = 8;
//  for(auto it = indexIterator(0,-1, VERTEX, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexface0[counter--],*it,"counter: " << counter);
//  }
//  WALBERLA_CHECK_EQUAL_2(counter, std::numeric_limits<size_t>::max());
//  counter = 8;
//  for(auto it = indexIterator(1,-1, VERTEX, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexface1[counter--],*it,"counter: " << counter);
//  }
//  WALBERLA_CHECK_EQUAL_2(counter, std::numeric_limits<size_t>::max());
//  counter = 8;
//  for(auto it = indexIterator(2,-1, VERTEX, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexface2[counter--],*it,"counter: " << counter);
//  }
//  WALBERLA_CHECK_EQUAL_2(counter, std::numeric_limits<size_t>::max());
//
//  /// CHECK GRAY CELL ITERATOR ///
//  std::vector<size_t> celledge0 = {45,46,47,48,49,50,51,52};
//  std::vector<size_t> celledge1 = {52,59,65,70,74,77,79,80};
//  std::vector<size_t> celledge2 = {80,78,75,71,66,60,53,45};
//  counter = 0;
//  for(auto it = indexIterator(0,1, CELL_GRAY, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(celledge0[counter++],*it,"counter: " << counter);
//  }
//  counter = 0;
//  for(auto it = indexIterator(1,1, CELL_GRAY, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(celledge1[counter++],*it,"counter: " << counter);
//  }
//  counter = 0;
//  for(auto it = indexIterator(2,1, CELL_GRAY, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(celledge2[counter++],*it,"counter: " << counter);
//  }
//  /// CHECK REVERSE ///
//  counter = 7;
//  for(auto it = indexIterator(0,-1, CELL_GRAY, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(celledge0[counter],*it,"counter: " << counter);
//    counter--;
//  }
//  counter = 7;
//  for(auto it = indexIterator(1,-1, CELL_GRAY, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(celledge1[counter],*it,"counter: " << counter);
//    counter--;
//  }
//  counter = 7;
//  for(auto it = indexIterator(2,-1, CELL_GRAY, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(celledge2[counter],*it,"counter: " << counter);
//    counter--;
//  }
//
//  /// CHECK BLUE CELL ITERATOR ///
//  std::vector<size_t> cellblueedge0 = {81,82,83,84,85,86,87};
//  std::vector<size_t> cellblueedge1 = {87,93,98,102,105,107,108};
//  std::vector<size_t> cellblueedge2 = {108,106,103,99,94,88,81};
//  counter = 0;
//  for(auto it = indexIterator(0,1, CELL_BLUE, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(cellblueedge0[counter],*it,"counter: " << counter);
//    counter++;
//  }
//  counter = 0;
//  for(auto it = indexIterator(1,1, CELL_BLUE, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(cellblueedge1[counter],*it,"counter: " << counter);
//    counter++;
//  }
//  counter = 0;
//  for(auto it = indexIterator(2,1, CELL_BLUE, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(cellblueedge2[counter],*it,"counter: " << counter);
//    counter++;
//  }
//  /// CHECK REVERSE ///
//  counter = 6;
//  for(auto it = indexIterator(0,-1, CELL_BLUE, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(cellblueedge0[counter],*it,"counter: " << counter);
//    counter--;
//  }
//  counter = 6;
//  for(auto it = indexIterator(1,-1, CELL_BLUE, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(cellblueedge1[counter],*it,"counter: " << counter);
//    counter--;
//  }
//  counter = 6;
//  for(auto it = indexIterator(2,-1, CELL_BLUE, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(cellblueedge2[counter],*it,"counter: " << counter);
//    counter--;
//  }
//  /// CHECK VERTEX_INNTER ITERATOR ///
//  std::vector<size_t> vertexinnerdge0 = {9,10,11,12,13,14,15,16};
//  std::vector<size_t> vertexinnerdge1 = {7,15,22,28,33,37,40,42};
//  std::vector<size_t> vertexinnerdge2 = {43,40,36,31,25,18,10,1};
//  counter = 0;
//  for(auto it = indexIterator(0,1, VERTEX_INNER, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexinnerdge0[counter],*it,"counter: " << counter);
//    counter++;
//  }
//  counter = 0;
//  for(auto it = indexIterator(1,1, VERTEX_INNER, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexinnerdge1[counter],*it,"counter: " << counter);
//    counter++;
//  }
//  counter = 0;
//  for(auto it = indexIterator(2,1, VERTEX_INNER, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexinnerdge2[counter],*it,"counter: " << counter);
//    counter++;
//  }
//  /// CHECK REVERSE ///
//  counter = 7;
//  for(auto it = indexIterator(0,-1, VERTEX_INNER, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexinnerdge0[counter],*it,"counter: " << counter);
//    counter--;
//  }
//  counter = 7;
//  for(auto it = indexIterator(1,-1, VERTEX_INNER, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexinnerdge1[counter],*it,"counter: " << counter);
//    counter--;
//  }
//  counter = 7;
//  for(auto it = indexIterator(2,-1, VERTEX_INNER, 3); it != indexIterator(); ++it){
//    WALBERLA_CHECK_EQUAL_3(vertexinnerdge2[counter],*it,"counter: " << counter);
//    counter--;
//  }
}
