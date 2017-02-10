#ifndef COMM_HPP
#define COMM_HPP

#include <mpi.h>

namespace hhg
{

class Comm
{
public:
  static Comm& get()
  {
    static Comm instance;

    MPI_Comm_size(MPI_COMM_WORLD, &instance.np);
    MPI_Comm_rank(MPI_COMM_WORLD, &instance.rk);

    return instance;
  }

  int rk;
  int np;

private:
  Comm() {}

  // Comm(Comm const&);
  // void operator=(Comm const&);

public:
  Comm(Comm const&) = delete;
  void operator=(Comm const&) = delete;
};

}

#endif /* COMM_HPP */