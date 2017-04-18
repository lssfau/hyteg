//
// Created by thoennes on 13.04.17.
//

#include <tinyhhg_core/tinyhhg.hpp>
#include <core/debug/TestSubsystem.h>
#include <core/mpi/MPIManager.h>

int main (int argc, char ** argv )
{

  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  hhg::Mesh mesh("./tri_1el.msh");

  return 0;
}
