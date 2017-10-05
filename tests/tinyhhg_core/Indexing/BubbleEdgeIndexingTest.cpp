#include "core/Environment.h"
#include "tinyhhg_core/tinyhhg.hpp"

//using hhg::P1BubbleFace::index;
using namespace hhg::BubbleEdge;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  /// CHECK VERTEX ///
  std::vector<size_t> refOne = {2,17,15,0,1,16};
  std::vector<size_t> refFive = {10,25,23,8,9,24};
  std::vector<size_t> result;
  for(auto n : neighbors)
  {
    size_t idx = indexFaceFromVertex<3>(1, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refOne.size(); ++i){
    WALBERLA_CHECK_EQUAL(refOne[i],result[i],"i: " << i);
  }
  result.clear();
  for(auto n : neighbors)
  {
    size_t idx = indexFaceFromVertex<3>(5, n);
    result.push_back(idx);
    //WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }
  for(size_t i = 0; i < refFive.size(); ++i){
    WALBERLA_CHECK_EQUAL(refFive[i],result[i],"i: " << i);
  }
  result.clear();
}
