#pragma once

#include "core/DataTypes.h"
#include "hyteg/types/flags.hpp"

namespace hyteg {

template < typename FunctionType >
class ProlongationOperator
{
public:
  virtual void prolongate( const FunctionType& function, const walberla::uint_t& sourceLevel, const DoFType& flag ) const = 0;

  virtual ~ProlongationOperator() = default;
};

} // namespace hyteg
