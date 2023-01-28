/*
* Copyright (c) 2023 Daniel Bauer.
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

#include <optional>
#include <string>
#include <vector>

#include "hyteg/eigen/typeAliases.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/types/PointND.hpp"

#include "../tests/hyteg/N1E1/common.hpp"
#include "KeyValueStore.hpp"

using namespace hyteg;

class Params
{
 public:
   std::string name;

   // system
   std::vector< real_t >                                               coefficients;
   n1e1::System                                                        system;
   std::optional< std::function< Eigen::Vector3r( const Point3D& ) > > initialGuess;

   // solver
   uint_t minLevel, maxLevel;
   bool   computeAndStoreLocalElementMatrices;
   // Chebyshev smoother
   uint_t chebyshevOrder;
   uint_t numSpectralRadiusEstIts;
   // hybrid smoother
   uint_t n1e1SmoothSteps, p1SmoothSteps;
   uint_t preSmoothSteps, postSmoothSteps;
   // stopping criteria
   uint_t                  nMaxIterations;
   std::optional< real_t > u2Reduction;
   std::optional< real_t > residual2Reduction;

   // output
   bool writeVTK;

   Params( const std::string _name )
   : name{ _name }
   , coefficients{ 1.0, 1.0 }
   , system{ n1e1::System::polynomialOnCube() }
   , initialGuess{}
   , minLevel{ 0 }
   , maxLevel{ 0 }
   , computeAndStoreLocalElementMatrices{ false }
   , chebyshevOrder{ 4 }
   , numSpectralRadiusEstIts{ 40 }
   , n1e1SmoothSteps{ 1 }
   , p1SmoothSteps{ 1 }
   , preSmoothSteps{ 1 }
   , postSmoothSteps{ 1 }
   , nMaxIterations{ 10 }
   , u2Reduction{}
   , residual2Reduction{}
   , writeVTK{ false }
   {}

   void store( KeyValueStore& store )
   {
      store.store( "/" + name + "/alpha", coefficients[0] );
      store.store( "/" + name + "/beta", coefficients[1] );
      store.store( "/" + name + "/minLevel", minLevel );
      store.store( "/" + name + "/maxLevel", maxLevel );
      store.store( "/" + name + "/computeAndStoreLocalElementMatrices", computeAndStoreLocalElementMatrices );
      store.store( "/" + name + "/chebyshevOrder", chebyshevOrder );
      store.store( "/" + name + "/numSpectralRadiusEstIts", numSpectralRadiusEstIts );
      store.store( "/" + name + "/n1e1SmoothSteps", n1e1SmoothSteps );
      store.store( "/" + name + "/p1SmoothSteps", p1SmoothSteps );
      store.store( "/" + name + "/preSmoothSteps", preSmoothSteps );
      store.store( "/" + name + "/postSmoothSteps", postSmoothSteps );
      store.store( "/" + name + "/nMaxIterations", nMaxIterations );

      if ( u2Reduction.has_value() )
      {
         store.store( "/" + name + "/u2Reduction", u2Reduction.value() );
      }
      if ( residual2Reduction.has_value() )
      {
         store.store( "/" + name + "/residual2Reduction", residual2Reduction.value() );
      }
   }
};

struct Results
{
   const uint_t numberOfGlobalDoFs;
   const real_t spectralRadius;
   const real_t initU2, finalU2;
   const real_t initResidual2, finalResidual2;
   const real_t finalErrL2;
   const uint_t nIterations;
};

Results solve( const Params& params );
Results smooth( const Params& params );
