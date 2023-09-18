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

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <unresolved_particles/data/shape/BaseShape.h>
#include <unresolved_particles/data/shape/Sphere.h>
#include <unresolved_particles/data/shape/HalfSpace.h>
#include <unresolved_particles/data/shape/CylindricalBoundary.h>
#include <unresolved_particles/data/shape/Box.h>
#include <unresolved_particles/data/shape/Ellipsoid.h>
#include <unresolved_particles/data/shape/ConvexPolyhedron.h>

#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

#include <memory>

namespace walberla {
namespace mpi {
   template< typename T,    // Element type of SendBuffer
             typename G>    // Growth policy of SendBuffer
   mpi::GenericSendBuffer<T,G>& operator<<( mpi::GenericSendBuffer<T,G> & buf,
                                            const std::shared_ptr<unresolved_particles::data::BaseShape>& bs )
   {
      buf.addDebugMarker( "up" );
      buf << bs->getShapeType();
      bs->pack(buf);
      return buf;
   }

   template< typename T>    // Element type  of RecvBuffer
   mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buf,
                                          std::shared_ptr<unresolved_particles::data::BaseShape>& bs )
   {
      using namespace unresolved_particles::data;

      buf.readDebugMarker( "up" );

      unresolved_particles::data::BaseShape::ShapeTypeT shapeType;
      buf >> shapeType;
      switch (shapeType)
      {
         case BaseShape::INVALID_SHAPE :
            bs = std::make_unique<unresolved_particles::data::BaseShape>();
            bs->unpack(buf);
            break;
         case Sphere::SHAPE_TYPE :
            bs = std::make_unique<unresolved_particles::data::Sphere>();
            bs->unpack(buf);
            break;
         case HalfSpace::SHAPE_TYPE :
            bs = std::make_unique<unresolved_particles::data::HalfSpace>();
            bs->unpack(buf);
            break;
         case CylindricalBoundary::SHAPE_TYPE :
            bs = std::make_unique<unresolved_particles::data::CylindricalBoundary>();
            bs->unpack(buf);
            break;
         case Box::SHAPE_TYPE :
            bs = std::make_unique<unresolved_particles::data::Box>();
            bs->unpack(buf);
            break;
         case Ellipsoid::SHAPE_TYPE :
            bs = std::make_unique<unresolved_particles::data::Ellipsoid>();
            bs->unpack(buf);
            break;
         case ConvexPolyhedron::SHAPE_TYPE :
            bs = std::make_unique<unresolved_particles::data::ConvexPolyhedron>();
            bs->unpack(buf);
            break;
         default : WALBERLA_ABORT("Shape type (" << shapeType << ") could not be determined!");
      }
      return buf;
   }
} //namespace mpi
} //namespace walberla