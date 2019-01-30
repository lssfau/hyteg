
#pragma once

#include "tinyhhg_core/gridtransferoperators/RestrictionOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"

namespace hhg {

class P2toP2QuadraticRestriction : public RestrictionOperator< P2Function< real_t > >
{
 public:
   inline void restrict( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      restrictAdditively( function, sourceLevel, flag );
   }

 private:
   void restrictWithPostCommunication( const P2Function< real_t >& function,
                                       const uint_t&               sourceLevel,
                                       const DoFType&              flag ) const;
   void restrictAdditively( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const;
};

} // namespace hhg