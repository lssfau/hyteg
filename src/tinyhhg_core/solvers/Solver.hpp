#pragma once

#include "core/DataTypes.h"

namespace hyteg {
template < class OperatorType >
class Solver
{
 public:
   /// solves the system A * x = b
   virtual void solve( const OperatorType&             A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const walberla::uint_t         level ) = 0;

   virtual ~Solver() = default;
};
} // namespace hyteg