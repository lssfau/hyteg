/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include <array>
#include <vector>

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/mesh/GmshReaderForMSH22.hpp"

namespace hyteg {

MeshInfo MeshInfo::fromGmshFile( const std::string& meshFileName )
{
   MeshInfo meshInfo;

   std::ifstream meshFile;
   meshFile.open( meshFileName.c_str() );

   WALBERLA_CHECK( !!meshFile, "[Mesh] Error opening file: " << meshFileName );

   std::string token;
   meshFile >> token; // $MeshFormat

   WALBERLA_CHECK_EQUAL( token, "$MeshFormat", "[Mesh] Missing: $MeshFormat" );

   meshFile >> token; // version

   meshFile.close();

   // delegate work to reader for specific MSH format
   if ( token == "2.2" )
   {
      meshInfo = GmshReaderForMSH22::readMesh( meshFileName );
   }
   else if ( token == "4.1" )
   {
      WALBERLA_ABORT( "Support for MSH 4.1 currently under construction" );
   }
   else
   {
      WALBERLA_ABORT( "MSH format " << token << " is not supported.\n"
                                    << "If you definitely need it, please open an issue on "
                                    << "https://i10git.cs.fau.de/hyteg/hyteg" );
   }

   return meshInfo;
}

} // namespace hyteg
