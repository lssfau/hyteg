/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <vector>

#include "core/Environment.h"
#include "core/timing/Timer.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/misc/dummy.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macrocell_vertexdof_to_vertexdof_add.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"

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

      std::map< hyteg::indexing::Index, double > stencil;
      for ( const auto & neighbor : hyteg::vertexdof::macrocell::neighborsWithCenter )
        stencil[hyteg::vertexdof::logicalIndexOffsetFromVertex( neighbor ) ] = walberla::real_c( std::rand() );

      double time(0.0);

      size_t iter = 2;
      while( timer.total() < 0.5 )
      {
         LIKWID_MARKER_START( "apply" );
         timer.reset();
         for( size_t i = 0; i < iter; ++i )
         {
            hyteg::vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add(
                dst.data(), src.data(), (int32_t) level, stencil );
            hyteg::misc::dummy( dst.data(), src.data() );
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