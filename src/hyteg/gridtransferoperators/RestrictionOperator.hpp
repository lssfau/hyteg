#pragma once

#include "core/DataTypes.h"
#include "hyteg/types/flags.hpp"

namespace hyteg {

template < typename FunctionType >
class RestrictionOperator
{
public:
  virtual void restrict( const FunctionType& function, const walberla::uint_t& sourceLevel, const DoFType& flag ) const = 0;

  virtual ~RestrictionOperator() = default;
};

} // namespace hyteg
