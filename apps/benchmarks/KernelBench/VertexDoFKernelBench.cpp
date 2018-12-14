#include <vector>
#include <iostream>

#include "core/timing/Timer.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/p1functionspace/generatedKernels/apply_2D_macroface_vertexdof_to_vertexdof_add.cpp"
//#include "tinyhhg_core/p1functionspace/generatedKernels/GeneratedKernels.hpp"

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT

   walberla::WcTimer timer;

   for( size_t level = 2; level < 14; ++level )
   {
      size_t edgeSize = (size_t) std::pow( 2u, level ) + 1;
      size_t faceSize = ( size_t )( edgeSize * ( edgeSize + 1u ) ) / 2;

      std::vector< double > src( faceSize );
      std::generate( src.begin(), src.end(), std::rand );
      std::vector< double > dst( faceSize );
      std::generate( dst.begin(), dst.end(), std::rand );
      std::vector< double > stencil( 7 );
      std::generate( stencil.begin(), stencil.end(), std::rand );

      timer.reset();

      hhg::vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add(
          dst.data(), src.data(), stencil.data(), (int64_t) level );

      timer.end();
      double time1 = timer.last();
      timer.reset();

      hhg::vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add_level_any(
          dst.data(), src.data(), stencil.data(), (int64_t) level );

      timer.end();
      double time2 = timer.last();

      std::cout << "special version: " << time1 << " any version: " << time2 << " speedup: " << time2 / time1 << std::endl;

   }
}