#include "tinyhhg_core/tinyhhg.hpp"

//using hhg::P1BubbleFace::index;
using namespace hhg::P1BubbleFace;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

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


  for(auto n : vertex_neighbors_with_center)
  {

    size_t idx = index<3>(6, 1, n);
    WALBERLA_LOG_INFO_ON_ROOT(enumStrings[n] << " " << idx);
  }

}
