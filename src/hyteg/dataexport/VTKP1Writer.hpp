/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#pragma once

#include "hyteg/dataexport/VTKOutput.hpp"

namespace hyteg {

class VTKOutput;

class VTKP1Writer
{
 public:
   static void write( const VTKOutput& mgr, std::ostream& output, const uint_t& level );

 private:
   template < typename value_t >
   static void writeScalarFunction( std::ostream&                                  output,
                                    const vertexdof::VertexDoFFunction< value_t >& function,
                                    const std::shared_ptr< PrimitiveStorage >&     storage,
                                    const uint_t&                                  level,
                                    bool                                           write2D,
                                    vtk::DataFormat                                vtkDataFormat );

   template < typename value_t >
   static void writeVectorFunction( std::ostream&                              output,
                                    const P1VectorFunction< value_t >&         function,
                                    const std::shared_ptr< PrimitiveStorage >& storage,
                                    const uint_t&                              level,
                                    bool                                       write2D,
                                    vtk::DataFormat                            vtkDataFormat );
};

} // namespace hyteg
