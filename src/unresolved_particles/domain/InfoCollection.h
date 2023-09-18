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

#include <blockforest/BlockForest.h>
#include <blockforest/loadbalancing/BlockInfo.h>
#include <blockforest/loadbalancing/InfoCollection.h>
#include <core/mpi/BufferSystem.h>

namespace walberla {
namespace unresolved_particles {
namespace domain {

template <typename Accessor>
void createWithNeighborhood(Accessor& ac, const BlockForest& bf, blockforest::InfoCollection& ic )
{
   using namespace walberla::blockforest;

   ic.clear();

   walberla::mpi::BufferSystem bs( MPIManager::instance()->comm(), 756 );

   for (size_t idx = 0; idx < ac.size(); ++idx)
   {
      for (auto blockIt = bf.begin(); blockIt != bf.end(); ++blockIt)
      {
         const blockforest::Block* block   = static_cast<const blockforest::Block*> (&(*blockIt));

         BlockInfo& info = ic[block->getId()];
         if (block->getAABB().contains(ac.getPosition(idx)))
         {
            if (data::particle_flags::isSet( ac.getFlags(idx), data::particle_flags::GHOST))
            {
               ++info.communicationWeight;
            } else
            {
               ++info.computationalWeight;
            }

            for (uint_t branchID = 0; branchID < 8; ++branchID)
            {
               const auto childID   = BlockID(block->getId(), branchID);
               const auto childAABB = bf.getAABBFromBlockId(childID);
               BlockInfo& childInfo = ic[childID];
               if (childAABB.contains(ac.getPosition(idx)))
               {
                  if (data::particle_flags::isSet( ac.getFlags(idx), data::particle_flags::GHOST))
                  {
                     ++childInfo.communicationWeight;
                  } else
                  {
                     ++childInfo.computationalWeight;
                  }
                  break; //particle can only be located within one child
               }
            }
            break; //particle can only be located within one block
         }
      }
   }

   for (auto blockIt = bf.begin(); blockIt != bf.end(); ++blockIt)
   {
      const blockforest::Block* block   = static_cast<const blockforest::Block*> (&(*blockIt));

      BlockInfo& info = ic[block->getId()];
      for( uint_t nb = uint_t(0); nb < block->getNeighborhoodSize(); ++nb )
      {
         bs.sendBuffer( block->getNeighborProcess(nb) ) << InfoCollection::value_type(block->getId(), info);
      }

      for (uint_t branchID = 0; branchID < 8; ++branchID)
      {
         const auto childID   = BlockID(block->getId(), branchID);

         BlockInfo& childInfo = ic[childID];

         for( uint_t nb = uint_t(0); nb < block->getNeighborhoodSize(); ++nb )
         {
            bs.sendBuffer( block->getNeighborProcess(nb) ) << InfoCollection::value_type(childID, childInfo);
         }
      }
   }

   // size of buffer is unknown and changes with each send
   bs.setReceiverInfoFromSendBufferState(false, true);
   bs.sendAll();

   for( auto recvIt = bs.begin(); recvIt != bs.end(); ++recvIt )
   {
      while( !recvIt.buffer().isEmpty() )
      {
         InfoCollectionPair val;
         recvIt.buffer() >> val;
         ic.insert(val);
      }
   }
}

} //namespace domain
} //namespace unresolved_particles
} //namespace walberla
