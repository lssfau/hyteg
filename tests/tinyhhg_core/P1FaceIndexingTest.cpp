#include "tinyhhg_core/p1functionspace/p1face.hpp"

#include "core/mpi/all.h"

//using hhg::P1BubbleFace::index;
using namespace hhg::P1Face;
//const Dir neighbors_with_center[] = {S, SE, W, C, E, NW, N};
int main(int argc, char* argv[])
{
  //this test is written for level 3
  walberla::mpi::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  /// CHECK VERTEX ///
  std::vector<size_t> refOneOne = {1,2,9,10,11,17,18};
  std::vector<size_t> refFiveTwo = {32,33,36,37,38,40,41};
  std::vector<size_t> result;
  for(auto n : neighbors_with_center)
  {
    size_t idx = index<3>(1, 1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refOneOne.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refOneOne[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : neighbors_with_center)
  {
    size_t idx = index<3>(5, 2, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refFiveTwo.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refFiveTwo[i],result[i],"i: " << i);
  }
  result.clear();

  /// CHECK VERTEX ITERATOR ///
  std::vector<size_t> vertexface0 = {0,1,2,3,4,5,6,7,8};
  std::vector<size_t> vertexface1 = {8,16,23,29,34,38,41,43,44};
  std::vector<size_t> vertexface2 = {44,42,39,35,30,24,17,9,0};
  uint_t counter = 0;
  for(auto it = indexIterator(0,1, VERTEX, 3); it != indexIterator(); ++it){
    WALBERLA_CHECK_EQUAL_3(vertexface0[counter++],*it,"counter: " << counter);
  }
  WALBERLA_CHECK_EQUAL_2(counter, 9);
  counter = 0;
  for(auto it = indexIterator(1,1, VERTEX, 3); it != indexIterator(); ++it){
    WALBERLA_CHECK_EQUAL_3(vertexface1[counter++],*it,"counter: " << counter);
  }
  WALBERLA_CHECK_EQUAL_2(counter, 9);
  counter = 0;
  for(auto it = indexIterator(2,1, VERTEX, 3); it != indexIterator(); ++it){
    WALBERLA_CHECK_EQUAL_3(vertexface2[counter++],*it,"counter: " << counter);
  }
  WALBERLA_CHECK_EQUAL_2(counter, 9);
  /// CHECK REVERSE ///
  counter = 8;
  for(auto it = indexIterator(0,-1, VERTEX, 3); it != indexIterator(); ++it){
    WALBERLA_CHECK_EQUAL_3(vertexface0[counter--],*it,"counter: " << counter);
  }
  WALBERLA_CHECK_EQUAL_2(counter, std::numeric_limits<size_t>::max());
  counter = 8;
  for(auto it = indexIterator(1,-1, VERTEX, 3); it != indexIterator(); ++it){
    WALBERLA_CHECK_EQUAL_3(vertexface1[counter--],*it,"counter: " << counter);
  }
  WALBERLA_CHECK_EQUAL_2(counter, std::numeric_limits<size_t>::max());
  counter = 8;
  for(auto it = indexIterator(2,-1, VERTEX, 3); it != indexIterator(); ++it){
    WALBERLA_CHECK_EQUAL_3(vertexface2[counter--],*it,"counter: " << counter);
  }
  WALBERLA_CHECK_EQUAL_2(counter, std::numeric_limits<size_t>::max());

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
