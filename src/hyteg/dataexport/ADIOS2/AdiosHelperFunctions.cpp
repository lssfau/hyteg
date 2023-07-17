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
#include "hyteg/dataexport/ADIOS2/AdiosHelperFunctions.hpp"

#include <adios2.h>

#include "core/logging/Logging.h"

#include "hyteg/dataexport/VTKOutput/VTKHelpers.hpp"

namespace hyteg::adiosHelpers {

// Note: Currently only the name "TIME" is allowed.
const std::string nameOfTimeStepVariable{ "TIME" };

std::string generateVTKMetaInfo( const std::vector< std::string >& namesOfPointDataFunctions,
                                 const std::vector< std::string >& namesOfCellDataFunctions )
{
   std::stringstream oStream;

   std::string myIndent1{ "  " };
   std::string myIndent2 = myIndent1 + myIndent1;
   std::string myIndent3 = myIndent1 + myIndent2;
   std::string myIndent4 = myIndent1 + myIndent3;

   // REMARK:
   //
   // Although we specify the names of the DataArrays in the vtk.xml part, these are
   // apparently not free to choose, but must be
   // - vertices
   // - connectivity
   // - types
   // oStream << R"(<?xml version="1.0"?>)" << '\n'
   oStream << R"(<VTKFile type="UnstructuredGrid" version="0.2" byte_order=")" << vtk::getByteOrder() << R"(">)" << '\n'
           << myIndent1 << R"(<UnstructuredGrid>)" << '\n'
           << myIndent2 << R"(<Piece NumberOfPoints="NumberOfVertices" NumberOfCells="NumberOfElements">)" << '\n'
           << myIndent3 << R"(<Points>)" << '\n'
           << myIndent4 << R"(<DataArray Name="vertices"/>)" << '\n'
           << myIndent3 << R"(</Points>)" << '\n'
           << myIndent3 << R"(<Cells>)" << '\n'
           << myIndent4 << R"(<DataArray Name="connectivity"/>)" << '\n'
           << myIndent4 << R"(<DataArray Name="types"/>)" << '\n'
           << myIndent3 << R"(</Cells>)" << '\n';

   // list of function names with node-based DoFs
   oStream << myIndent3 << R"(<PointData>)" << '\n';
   for ( const auto& name : namesOfPointDataFunctions )
   {
      oStream << myIndent4 << R"(<DataArray Name=")" << name << R"("/>)" << '\n';
   }
   // seems we always need this?
   oStream << myIndent4 << R"(<DataArray Name="TIME">)" << nameOfTimeStepVariable << R"(</DataArray>)" << '\n';
   oStream << myIndent3 << R"(</PointData>)" << '\n';

   // REMARK:
   //
   // cell data functions seem to currently be ignored by Paraview's BP reader

   // list of function names with cell-based DoFs
   oStream << myIndent3 << R"(<CellData>)" << '\n';
   for ( const auto& name : namesOfCellDataFunctions )
   {
      oStream << myIndent4 << R"(<DataArray Name=")" << name << R"("/>)" << '\n';
   }
   oStream << myIndent3 << R"(</CellData>)" << '\n';

   // wrap-up
   oStream << myIndent2 << R"(</Piece>)" << '\n' << myIndent1 << R"(</UnstructuredGrid>)" << '\n' << R"(</VTKFile>)";

   return oStream.str();
};

} // namespace hyteg::adiosHelpers
