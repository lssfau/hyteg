#pragma once

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/gridtransferoperators/RestrictionOperator.hpp"

namespace hyteg {

class P2P1StokesToP2P1StokesRestriction : public RestrictionOperator< P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2toP2QuadraticRestriction VelocityRestriction_T;
   typedef P1toP1LinearRestriction    PressureRestriction_T;

   P2P1StokesToP2P1StokesRestriction()
   : projectMeanAfterRestriction_( false )
   {}
   P2P1StokesToP2P1StokesRestriction( bool projectMeanAfterRestriction )
   : projectMeanAfterRestriction_( projectMeanAfterRestriction )
   {}

   void
       restrict( const P2P1TaylorHoodFunction< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
   {
      quadraticRestrictionOperator_.restrict( function.u, sourceLevel, flag );
      quadraticRestrictionOperator_.restrict( function.v, sourceLevel, flag );
      quadraticRestrictionOperator_.restrict( function.w, sourceLevel, flag );
      linearRestrictionOperator_.restrict( function.p, sourceLevel, flag );

      if ( projectMeanAfterRestriction_ )
      {
         vertexdof::projectMean( function.p, sourceLevel - 1 );
      }
   }

 private:
   P2toP2QuadraticRestriction quadraticRestrictionOperator_;
   P1toP1LinearRestriction    linearRestrictionOperator_;

   bool projectMeanAfterRestriction_;
};
} // namespace hyteg