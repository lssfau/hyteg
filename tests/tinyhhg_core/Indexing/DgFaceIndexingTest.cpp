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

  /// CHECK VERTEX ///
  std::vector <size_t> refOneOne = {1, 9, 8, 37, 43, 36};
  std::vector <size_t> refTwoFive = {28, 32, 31, 60, 62, 59};
  std::vector <size_t> result;
}