
#pragma once

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"

namespace hyteg {

class P1P1StokesToP1P1StokesRestriction : public RestrictionOperator< P1StokesFunction< real_t > >
{
 public:
   P1P1StokesToP1P1StokesRestriction()
   : projectMeanAfterRestriction_( false )
   {}
   P1P1StokesToP1P1StokesRestriction( bool projectMeanAfterRestriction )
   : projectMeanAfterRestriction_( projectMeanAfterRestriction )
   {}

   typedef P1toP1LinearRestriction VelocityRestriction_T;
   typedef P1toP1LinearRestriction PressureRestriction_T;

   void restrict( const P1StokesFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      restrictionOperator_.restrict( function.u, sourceLevel, flag );
      restrictionOperator_.restrict( function.v, sourceLevel, flag );
      restrictionOperator_.restrict( function.w, sourceLevel, flag );
      restrictionOperator_.restrict( function.p, sourceLevel, flag );

      if ( projectMeanAfterRestriction_ )
      {
         vertexdof::projectMean( function.p, sourceLevel - 1 );
      }
   }

 private:
   P1toP1LinearRestriction restrictionOperator_;

   bool projectMeanAfterRestriction_;
};
} // namespace hyteg