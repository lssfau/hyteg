#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "mesh.hpp"
#include <core/all.h>

namespace hhg
{

class Operator
{
public:
  Operator(Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
    : mesh(_mesh), minLevel(_minLevel), maxLevel(_maxLevel), memory_id(std::numeric_limits<std::size_t>::max()), rank(walberla::uint_c(walberla::mpi::MPIManager::instance()->rank()) )
  {
  }

  virtual ~Operator()
  {
  }

  Mesh& mesh;
  size_t minLevel;
  size_t maxLevel;
  size_t memory_id;
  size_t rank;
};

}

#endif /* OPERATOR_HPP */
