
#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"

namespace hhg {

class P1toP1LinearProlongation
{
 public:
   inline void operator()( const P1Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const
   {
      if ( function.isDummy() )
        return;

      if ( function.getStorage()->hasGlobalCells() )
      {
        prolongate3D( function, sourceLevel, flag );
      }
      else
      {
        prolongate2D( function, sourceLevel, flag );
      }
   }

 private:

   void prolongate2D( const P1Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const;

   void prolongate3D( const P1Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const;

   void prolongateMacroVertex2D( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const;

   void prolongateMacroEdge2D( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const;

   void prolongateMacroFace2D( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const;

};

} // namespace hhg