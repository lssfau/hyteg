#include "tinyhhg_core/p1functionspace/P1Face.hpp"

#include "core/mpi/all.h"

using namespace hhg::P1Face;

void checkIndices(uint_t col, uint_t row, std::vector<uint_t> ref, uint_t type){
  std::vector<size_t> result;
  switch(type){
    //vertex
    case 0:
      for(auto n : CoordsVertex::neighbors_with_center)
      {
        result.push_back(CoordsVertex::index<3>(col, row, n));
      }
      break;
    case 1:
      for(auto n : CoordsCellGray::neighbors)
      {
        result.push_back(CoordsCellGray::index<3>(col, row, n));
      }
      break;
    case 2:
      for(auto n : CoordsCellBlue::neighbors)
      {
        result.push_back(CoordsCellBlue::index<3>(col, row, n));
      }
      break;
    default:
      WALBERLA_ABORT("wrong type")
  }

  for(size_t i = 0; i < ref.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(ref[i],result[i],"col: " << col << " row: " << row <<" i: " << i);
  }
}

//const Dir neighbors_with_center[] = {S, SE, W, C, E, NW, N};
int main(int argc, char* argv[])
{
  //this test is written for level 3
  walberla::mpi::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  /// CHECK VERTEX INDEXING ///
  std::vector<size_t> refOneOne = {10,1,2,11,18,17,9};
  checkIndices(1,1,refOneOne,0);
  std::vector<size_t> refFiveTwo = {22,14,15,23,29,28,21};
  checkIndices(5,2,refFiveTwo,0);

  /// CHECK CELL GRAY INDEXING ///
  checkIndices(2,2,{19,20,26},1);
  checkIndices(2,4,{32,33,37},1);

  /// CHECK CELL BLUE INDEXING ///
  checkIndices(3,3,{28,33,34},2);
  checkIndices(4,1,{14,21,22},2);


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
