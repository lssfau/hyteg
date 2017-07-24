#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "mesh.hpp"
#include <core/all.h>

namespace hhg
{

class Operator
{
public:
  Operator(PrimitiveStorage& storage, size_t _minLevel, size_t _maxLevel)
    : storage_(storage), minLevel(_minLevel), maxLevel(_maxLevel), memory_id(std::numeric_limits<std::size_t>::max())
  {
  }

  virtual ~Operator()
  {
  }

  const PrimitiveStorage& storage_;
  const uint_t minLevel;
  const uint_t maxLevel;
  size_t memory_id;
};

}

#endif /* OPERATOR_HPP */
