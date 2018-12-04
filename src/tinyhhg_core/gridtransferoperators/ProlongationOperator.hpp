#pragma once

#include "core/DataTypes.h"
#include "tinyhhg_core/types/flags.hpp"

namespace hhg {

template < typename FunctionType >
class ProlongationOperator
{
public:
  virtual void prolongate( const FunctionType& function, const walberla::uint_t& sourceLevel, const DoFType& flag ) = 0;

  virtual ~ProlongationOperator() = default;
};

} // namespace hhg
