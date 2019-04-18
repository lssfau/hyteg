#include <iostream>
#include <vector>

#include "core/Environment.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/misc/dummy.hpp"
#include "tinyhhg_core/p1functionspace/generatedKernels/apply_3D_macrocell_vertexdof_to_vertexdof_add.cpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"

int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;
   LIKWID_MARKER_REGISTER( "apply" );

   walberla::WcTimer timer;

   for( size_t level = 9; level < 10; ++level )
   {
      size_t edgeSize = (size_t) std::pow( 2u, level ) + 1;
      //size_t faceSize = ( size_t )( edgeSize * ( edgeSize + 1u ) ) / 2;
      size_t tetSize  = ( size_t )( ( edgeSize + 2 ) * ( edgeSize + 1 ) * edgeSize ) / 6;

      std::vector< double > src( tetSize );
      std::generate( src.begin(), src.end(), std::rand );
      std::vector< double > dst( tetSize );
      std::generate( dst.begin(), dst.end(), std::rand );

      hhg::vertexdof::macrocell::StencilMap_T stencil;
      for ( const auto & neighbor : hhg::vertexdof::macrocell::neighborsWithCenter )
        stencil[ hhg::vertexdof::logicalIndexOffsetFromVertex( neighbor ) ] = walberla::real_c( std::rand() );

      double time(0.0);

      size_t iter = 2;
      while( timer.total() < 0.5 )
      {
         LIKWID_MARKER_START( "apply" );
         timer.reset();
         for( size_t i = 0; i < iter; ++i )
         {

            hhg::vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add(
                dst.data(), src.data(), (int64_t) level, stencil );
            hhg::misc::dummy( dst.data(), src.data() );
         }
         timer.end();
         LIKWID_MARKER_STOP( "apply" );


         time = timer.total();
         iter *= 2;
      }
      iter /= 2;

      timer.reset();


      std::cout << "Level: " << level << " time per iteration: " << time / (double) iter << std::endl;

      LIKWID_MARKER_CLOSE;
   }
}