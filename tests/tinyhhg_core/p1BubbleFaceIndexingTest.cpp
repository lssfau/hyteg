#include "tinyhhg_core/tinyhhg.hpp"

//using hhg::P1BubbleFace::index;
using namespace hhg::P1BubbleFace;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  std::vector<real_t> refOneOne = {10,1,2,11,18,17,9,46,54,53,82,88,81};
  std::string enumStrings[] = {
      "VERTEX_C",
      "VERTEX_S",
      "VERTEX_SE",
      "VERTEX_E",
      "VERTEX_N",
      "VERTEX_NW",
      "VERTEX_W",
      "CELL_GRAY_SE",
      "CELL_GRAY_NE",
      "CELL_GRAY_NW",
      "CELL_BLUE_SE",
      "CELL_BLUE_NW",
      "CELL_BLUE_SW"
  };

  std::vector<real_t> result;
  for(auto n : vertex_neighbors_with_center)
  {

    size_t idx = indexVertex<3>(1, 1, n);
    result.push_back(walberla::real_c(idx));
    WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }

  for(size_t i = 0; i < refOneOne.size(); ++i){
    WALBERLA_CHECK_EQUAL_3(refOneOne[i],result[i],"i: " << i);
  }
}
