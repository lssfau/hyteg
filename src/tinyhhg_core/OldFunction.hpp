#pragma once

#include "mesh.hpp"
#include "operator.hpp"
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"

#include <string>
#include <functional>

namespace hhg
{

class OldFunction
{
public:

  OldFunction(const std::string& _name, Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
    : name(_name), mesh(_mesh), minLevel(_minLevel), maxLevel(_maxLevel), memory_id(std::numeric_limits<std::size_t>::max()), rank(walberla::uint_c(walberla::mpi::MPIManager::instance()->rank()))
  {
  }

  virtual ~OldFunction()
  {
  }

  std::string name;
  Mesh& mesh;
  size_t minLevel;
  size_t maxLevel;
  size_t memory_id;
  size_t rank;
};

}
