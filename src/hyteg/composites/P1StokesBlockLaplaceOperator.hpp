#pragma once

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

class P1StokesBlockLaplaceOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   P1StokesBlockLaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   {}

   void apply( const P1StokesFunction< real_t >& src,
               const P1StokesFunction< real_t >& dst,
               const size_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType ) const
   {
      A.apply( src.u, dst.u, level, flag, updateType );
      A.apply( src.v, dst.v, level, flag, updateType );
      A.apply( src.p, dst.p, level, flag, updateType );
   }

   P1ConstantLaplaceOperator A;
};

} // namespace hyteg
