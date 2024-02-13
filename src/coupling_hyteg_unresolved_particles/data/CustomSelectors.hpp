/*
* Copyright (c) 2023 Nils Kohl.
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

#include "unresolved_particles/vtk/OutputSelector.h"

namespace walberla {
namespace unresolved_particles {
namespace data {

class SelectParticleCustomRealElement
{
 public:
   SelectParticleCustomRealElement( uint_t element )
   : element_( element )
   {}

   using return_type = walberla::real_t;

   walberla::real_t& operator()( data::Particle& p ) const
   {
      WALBERLA_ASSERT_GREATER( p.getCustomRealRef().size(), element_ );
      return p.getCustomRealRef()[element_];
   }

   walberla::real_t& operator()( data::Particle&& p ) const
   {
      WALBERLA_ASSERT_GREATER( p.getCustomRealRef().size(), element_ );
      return p.getCustomRealRef()[element_];
   }

   walberla::real_t const& operator()( const data::Particle& p ) const
   {
      WALBERLA_ASSERT_GREATER( p.getCustomReal().size(), element_ );
      return p.getCustomReal()[element_];
   }

   uint_t element_;
};

class SelectParticleCustomIntElement
{
 public:
   SelectParticleCustomIntElement( uint_t element )
   : element_( element )
   {}

   using return_type = int;

   int& operator()( data::Particle& p ) const
   {
      WALBERLA_ASSERT_GREATER( p.getCustomIntRef().size(), element_ );
      return p.getCustomIntRef()[element_];
   }

   int& operator()( data::Particle&& p ) const
   {
      WALBERLA_ASSERT_GREATER( p.getCustomIntRef().size(), element_ );
      return p.getCustomIntRef()[element_];
   }

   int const& operator()( const data::Particle& p ) const
   {
      WALBERLA_ASSERT_GREATER( p.getCustomInt().size(), element_ );
      return p.getCustomInt()[element_];
   }

   uint_t element_;
};

} // namespace data
} // namespace unresolved_particles
} // namespace walberla