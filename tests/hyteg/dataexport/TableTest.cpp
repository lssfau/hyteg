/*
 * Copyright (c) 2023 Marcus Mohr.
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

#include "hyteg/dataexport/LaTeX/Table.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export for 2D meshes:" );

   const uint_t nSamples = 10;

   hyteg::latex::Table< 3 > table( { R"(\varphi)", R"(sin(\varphi))", R"(cos(\varphi))" } );

   for ( uint_t k = 0; k < nSamples / 2; ++k )
   {
      real_t phi = real_c( k ) * pi / real_c( nSamples - 1 );
      table.addElement( k, 0, phi * real_c( 180 ) / pi );
      table.addElement( k, 1, std::sin( phi ) );
      table.addElement( k, 2, std::cos( phi ) );
   }

   for ( uint_t k = nSamples / 2; k < nSamples; ++k )
   {
      real_t phi = real_c( k ) * pi / real_c( nSamples - 1 );
      table.pushRow( phi * real_c( 180 ) / pi, std::sin( phi ), std::cos( phi ) );
   }

   table.write( ".", "TableTest" );
}
