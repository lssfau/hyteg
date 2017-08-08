#include "tinyhhg_core/p1functionspace/P1Face.hpp"

#include "core/mpi/all.h"

using namespace hhg::P1Face;

void checkIndices(uint_t col, uint_t row, std::vector<uint_t> ref){
  std::vector<size_t> result;
  for(auto n : neighbors_with_center)
  {
    size_t idx = index<3>(col, row, n);
    result.push_back(idx);
  }
  for(size_t i = 0; i < ref.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(ref[i],result[i],"i: " << i);
  }
}

//const Dir neighbors_with_center[] = {S, SE, W, C, E, NW, N};
int main(int argc, char* argv[])
{
  //this test is written for level 3
  walberla::mpi::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  /// CHECK VERTEX ///
  std::vector<size_t> refOneOne = {1,2,9,10,11,17,18};
  checkIndices(1,1,refOneOne);
  std::vector<size_t> refFiveTwo = {14,15,21,22,23,28,29};
  checkIndices(5,2,refFiveTwo);

  /// CHECK VERTEX ITERATOR ///

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
