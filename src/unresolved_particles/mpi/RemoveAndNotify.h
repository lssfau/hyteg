//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <unresolved_particles/data/ParticleStorage.h>
#include <core/mpi/BufferSystem.h>

namespace walberla {
namespace unresolved_particles {
namespace mpi {

data::ParticleStorage::iterator removeAndNotify( walberla::mpi::BufferSystem& bs,
                                                 data::ParticleStorage& ps,
                                                 data::ParticleStorage::iterator& pIt );

}  // namespace mpi
}  // namespace unresolved_particles
}  // namespace walberla
