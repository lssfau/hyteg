/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "core/singleton/Singleton.h"
#include "core/OpenMP.h"

namespace walberla {

/// \brief Class that may be useful to globally manage OpenMP functionality.
class OpenMPManager : public walberla::singleton::Singleton< OpenMPManager >
{
 public:

   WALBERLA_BEFRIEND_SINGLETON;

   /// \brief Forces the number of OpenMP threads to 1 after this call.
   void forceSerial()
   {
      omp_set_num_threads(1);
   }

   /// \brief Resets the (maximum) number of OpenMP threads to the number of threads that OpenMP was configured to
   ///        during construction of this manager.
   void resetToParallel()
   {
      omp_set_num_threads(maxNumThreads_);
   }

   /// \brief Returns the current (maximum) number of threads.
   int numThreads() const
   {
      return omp_get_max_threads();
   }

 private:

   OpenMPManager() : maxNumThreads_( omp_get_max_threads() ) {}

   int maxNumThreads_;

};

}