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
#include "constantStencilOperator/P1generatedKernels/apply_2D_macroface_vertexdof_to_vertexdof_add.cpp"

int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;
   LIKWID_MARKER_REGISTER( "apply_special" );
   LIKWID_MARKER_REGISTER( "apply_any" );

   walberla::WcTimer timerSpecial;
   walberla::WcTimer timerAny;

   for( size_t level = 14; level < 15; ++level )
   {
      size_t edgeSize = (size_t) std::pow( 2u, level ) + 1;
      size_t faceSize = ( size_t )( edgeSize * ( edgeSize + 1u ) ) / 2;

      std::vector< double > src( faceSize );
      std::generate( src.begin(), src.end(), std::rand );
      std::vector< double > dst( faceSize );
      std::generate( dst.begin(), dst.end(), std::rand );
      std::vector< double > stencil( 7 );
      std::generate( stencil.begin(), stencil.end(), std::rand );

      double timeSpecial(0.0);
      double timeAny(0.0);

      size_t iter = 2;
      while( timerSpecial.total() < 0.5 )
      {
         LIKWID_MARKER_START( "apply_special" );
         timerSpecial.reset();
         for( size_t i = 0; i < iter; ++i )
         {
            hyteg::vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add(
                dst.data(), src.data(), stencil.data(), (int32_t) level );
            hyteg::misc::dummy( dst.data(), src.data() );
         }
         timerSpecial.end();
         LIKWID_MARKER_STOP( "apply_special" );

         LIKWID_MARKER_START( "apply_any" );
         timerAny.reset();
         for( size_t i = 0; i < iter; ++i )
         {
            hyteg::vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add_level_any(
                dst.data(), src.data(), stencil.data(), (int32_t) level );
            hyteg::misc::dummy( dst.data(), src.data() );
         }
         timerAny.end();
         LIKWID_MARKER_STOP( "apply_any" );

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

      LIKWID_MARKER_CLOSE;
   }
}