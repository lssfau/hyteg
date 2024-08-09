/*
 * Copyright (c) 2024 Benjamin Mann, Marcus Mohr.
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

#include "core/mpi/Environment.h"

#include "hyteg/adaptiverefinement/mesh.hpp"

namespace hyteg {

/**
 * \page BA.05_AMR Tutorial BA.05 - Applying adaptive mesh refinement (red-green)
 *
 * \dontinclude tutorials/basics-for-apps/BA.05_AMR.cpp
 *
 * \brief In this tutorial we will set up a mesh that can be adaptively refined using the class adaptiveRefinement::Mesh.
 *
 * \section T05-intro Introduction
 *
 * In order to apply adaptive refinement, we require an additional data structure between SetupPrimitiveStorage and PrimitiveStorage, containing the necessary tree structure with the entire refinement history.
 *
 * The process introduced in \ref BA.01_PrimitiveStorage is therefore extended by another step, yielding four steps in total:
 *
 * -# creating a MeshInfo object from a mesh file
 * -# creating a SetupPrimitiveStorage from the MeshInfo object
 * -# creating an adaptiveRefinement::Mesh from the SetupPrimitiveStorage object
 * -# creating a distributed PrimitiveStorage from the current state of the refined mesh
 *
 * \section BA05-setupstorage Creating MeshInfo and SetupPrimitiveStorage
 *
 * Creating the initial mesh requires a SetupPrimitiveStorage, that is created as introduced in \ref BA.01_PrimitiveStorage.
 *
 * \snippet{trimleft} this SetupPrimitiveStorage
 *
 * \section BA05-adaptivemesh Creating an adaptive mesh
 *
 * We now use the SetupPrimitiveStorage to create the initial grid, instead of directly creating a PrimitiveStorage.
 *
 * \snippet{trimleft} this AdaptiveMesh
 *
 * As with the regular approach from \ref BA.01_PrimitiveStorage, we require load balancing.
 * Instead of calling this on the Setupstorage, we can directly call a loadbalancer from the adaptive mesh.
 * Again, there are some simple load balancers available in the library, defaulting to round robin.
 *
 * \snippet{trimleft} this Loadbalancing
 *
 * Finally we are able to create the distributed PrimitiveStorage instance.
 * We do this by calling the method make_storage()
 *
 * \snippet{trimleft} this PrimitiveStorage
 *
 * \section BA05-refine Refining the mesh
 *
 * First, we need to select the elements that shall be refined.
 * Note that, since we apply red-green refinement, in general, more elements than the ones you selected will be refined.
 * Here, we simply mark the first macro element for refinement.
 * There are more sophisticated marking strategies available.
 * These will be introduced in a separate tutorial about error estimators.
 * Note that, after applying refinement, we must apply loadbalancing again.
 *
 * \snippet{trimleft} this RefineRG
 *
 * Once the mesh is refined, we can create a new PrimitiveStorage corresponding to the refined mesh and work on this new storage as usual.
 *
 * \snippet{trimleft} this RefinedPrimitiveStorage
 *
 * \section BA05-Code Complete Program
 * \include tutorials/basics-for-apps/BA.05_AMR.cpp
 *
 *
 */

void amrTutorial()
{
   /// [SetupPrimitiveStorage]
   auto                  meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/tri_2el.msh" ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   /// [SetupPrimitiveStorage]

   /// [AdaptiveMesh]
   adaptiveRefinement::Mesh adaptive_mesh( setupStorage );
   /// [AdaptiveMesh]

   /// [Loadbalancing]
   adaptive_mesh.loadbalancing();
   /// [Loadbalancing]

   /// [PrimitiveStorage]
   std::shared_ptr< hyteg::PrimitiveStorage > storage = adaptive_mesh.make_storage();
   /// [PrimitiveStorage]

   /// [RefineRG]
   // mark one element for refinement
   auto toRefine = setupStorage.getFaces().begin()->first;
   adaptive_mesh.refineRG( { toRefine } );
   // apply load balancing for the resulting mesh
   adaptive_mesh.loadbalancing();
   /// [RefineRG]

   /// [RefinedPrimitiveStorage]
   storage = adaptive_mesh.make_storage();
   /// [RefinedPrimitiveStorage]
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   hyteg::amrTutorial();
   return 0;
}
