#pragma once

#include "core/DataTypes.h"

namespace hhg {
template < class OperatorType >
class Solver
{
 public:
   /// solves the system A * x = b
   virtual void solve( const OperatorType&             A,
                       typename OperatorType::srcType& x,
                       typename OperatorType::dstType& b,
                       const walberla::uint_t&         level ) = 0;

   virtual ~Solver() = default;
};
} // namespace hhg