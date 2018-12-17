#include <iostream>
#include <vector>

#include "core/timing/Timer.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/misc/dummy.hpp"
#include "tinyhhg_core/p1functionspace/generatedKernels/apply_2D_macroface_vertexdof_to_vertexdof_add.cpp"
//#include "tinyhhg_core/p1functionspace/generatedKernels/GeneratedKernels.hpp"

int main( int, char** )
{
   walberla::WcTimer timerSpecial;
   walberla::WcTimer timerAny;

   for( size_t level = 2; level < 15; ++level )
   {
      size_t edgeSize = (size_t) std::pow( 2u, level ) + 1;
      size_t faceSize = ( size_t )( edgeSize * ( edgeSize + 1u ) ) / 2;

      std::vector< double > src( faceSize );
      std::generate( src.begin(), src.end(), std::rand );
      std::vector< double > dst( faceSize );
      std::generate( dst.begin(), dst.end(), std::rand );
      std::vector< double > stencil( 7 );
      std::generate( stencil.begin(), stencil.end(), std::rand );

      double timeSpecial, timeAny;

      size_t iter = 2;
      while( timerSpecial.total() < 0.5 )
      {
         timerSpecial.reset();
         for( size_t i = 0; i < iter; ++i )
         {
            hhg::vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add(
                dst.data(), src.data(), stencil.data(), (int64_t) level );
            hhg::misc::dummy( dst.data(), src.data() );
         }
         timerSpecial.end();

         timerAny.reset();
         for( size_t i = 0; i < iter; ++i )
         {
            hhg::vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add_level_any(
                dst.data(), src.data(), stencil.data(), (int64_t) level );
            hhg::misc::dummy( dst.data(), src.data() );
         }
         timerAny.end();

         timeSpecial = timerSpecial.total();
         timeAny     = timerAny.total();
         iter *= 2;
      }
      iter /= 2;

      timerSpecial.reset();
      timerAny.reset();

      std::cout << "Level: " << level << " special version: " << timeSpecial / (double) iter
                << " any version: " << timeAny / (double) iter << " speedup (Any/special): " << timeAny / timeSpecial
                << std::endl;
   }
}