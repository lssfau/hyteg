
#pragma once

#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

class P2toP2QuadraticRestriction : public RestrictionOperator< P2Function< real_t > >
{
 public:
   inline void restrict( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      if ( function.isDummy() )
         return;

      if ( function.getStorage()->hasGlobalCells() )
      {
         restrictAdditively3D( function, sourceLevel, flag );
      }
      else
      {
         restrictAdditively( function, sourceLevel, flag );
      }
   }

 private:
   void restrictWithPostCommunication( const P2Function< real_t >& function,
                                       const uint_t&               sourceLevel,
                                       const DoFType&              flag ) const;
   void restrictAdditively( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const;
   void restrictAdditively3D( const P2Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const;
};

} // namespace hyteg