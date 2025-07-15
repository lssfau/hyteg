/*
* Copyright (c) 2025 Andreas Burkhart
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
#include "ElementwiseDoFValueOperator.hpp"

namespace hyteg {

void communicateDoFValue( const P1Function< real_t >& f, const uint_t level )
{
   if ( f.getStorage()->hasGlobalCells() )
   {
      f.communicate< Face, Cell >( level );
      f.communicate< Edge, Cell >( level );
      f.communicate< Vertex, Cell >( level );
   }
   else
   {
      communication::syncFunctionBetweenPrimitives( f, level );
   }
}

void communicateDoFValue( const P2Function< real_t >& f, const uint_t level )
{
   communicateDoFValue( f.getVertexDoFFunction(), level );
   communicateDoFValue( f.getEdgeDoFFunction(), level );
}

void communicateDoFValue( const EdgeDoFFunction< real_t >& f, const uint_t level )
{
   if ( f.getStorage()->hasGlobalCells() )
   {
      f.communicate< Face, Cell >( level );
      f.communicate< Edge, Cell >( level );
      f.communicate< Vertex, Cell >( level );
   }
   else
   {
      communication::syncFunctionBetweenPrimitives( f, level );
   }
}

} // namespace hyteg