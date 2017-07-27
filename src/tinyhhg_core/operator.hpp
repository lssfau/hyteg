#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "mesh.hpp"
#include <core/all.h>

namespace hhg
{

class Operator
{
public:
  Operator(const std::shared_ptr<PrimitiveStorage> & storage, uint_t minLevel, uint_t maxLevel)
    : storage_(storage), minLevel_(minLevel), maxLevel_(maxLevel)
  {
  }

  virtual ~Operator()
  {
  }

 protected:
  const std::shared_ptr< PrimitiveStorage > storage_;
  const uint_t minLevel_;
  const uint_t maxLevel_;
};

}

#endif /* OPERATOR_HPP */
