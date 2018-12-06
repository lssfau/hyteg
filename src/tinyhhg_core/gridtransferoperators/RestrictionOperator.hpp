#pragma once

#include "core/DataTypes.h"
#include "tinyhhg_core/types/flags.hpp"

namespace hhg {

template < typename FunctionType >
class RestrictionOperator
{
public:
  virtual void restrict( const FunctionType& function, const walberla::uint_t& sourceLevel, const DoFType& flag ) const = 0;

  virtual ~RestrictionOperator() = default;
};

} // namespace hhg
