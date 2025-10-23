/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include <fstream>

#include "core/DataTypes.h"
#include "core/config/Config.h"

#include "../Parameters/NondimensionalisationParameters.hpp"

using walberla::real_t;

namespace MantleConvection {

// expects rangeX to be sorted
real_t linearInterpolation( real_t x, const std::vector< real_t >& rangeX, const std::vector< real_t >& values )
{
   real_t factorX = ( x - rangeX.front() ) / ( rangeX.back() - rangeX.front() );
   real_t factorXInner;
   uint_t indexX0;
   uint_t indexX1;
   if ( factorX >= real_c( 1.0 ) )
   {
      indexX0      = rangeX.size() - 1;
      indexX1      = indexX0;
      factorXInner = real_c( 1.0 );
   }
   else if ( factorX <= real_c( 0.0 ) )
   {
      indexX0      = 0;
      indexX1      = indexX0;
      factorXInner = real_c( 0.0 );
   }
   else
   {
      auto iX      = std::lower_bound( rangeX.begin(), rangeX.end(), x );
      indexX1      = static_cast< uint_t >( std::distance( rangeX.begin(), iX ) );
      indexX0      = indexX1 - 1;
      factorXInner = ( x - rangeX.at( indexX0 ) ) / ( rangeX.at( indexX1 ) - rangeX.at( indexX0 ) );
   }

   // linear interpolation in x direction
   return values.at( indexX0 ) + factorXInner * ( values.at( indexX1 ) - values.at( indexX0 ) );
}

// expects rangeX, rangeY to be sorted
real_t bilinearInterpolation( real_t                                      x,
                              real_t                                      y,
                              const std::vector< real_t >&                rangeX,
                              const std::vector< real_t >&                rangeY,
                              const std::vector< std::vector< real_t > >& values )
{
   // x direction
   real_t factorX = ( x - rangeX.front() ) / ( rangeX.back() - rangeX.front() );
   real_t factorXInner;
   uint_t indexX0;
   uint_t indexX1;
   if ( factorX >= real_c( 1.0 ) )
   {
      indexX0      = rangeX.size() - 1;
      indexX1      = indexX0;
      factorXInner = real_c( 1.0 );
   }
   else if ( factorX <= real_c( 0.0 ) )
   {
      indexX0      = 0;
      indexX1      = indexX0;
      factorXInner = real_c( 0.0 );
   }
   else
   {
      auto iX      = std::lower_bound( rangeX.begin(), rangeX.end(), x );
      indexX1      = static_cast< uint_t >( std::distance( rangeX.begin(), iX ) );
      indexX0      = indexX1 - 1;
      factorXInner = ( x - rangeX.at( indexX0 ) ) / ( rangeX.at( indexX1 ) - rangeX.at( indexX0 ) );
   }

   // y direction
   real_t factorY = ( y - rangeY.front() ) / ( rangeY.back() - rangeY.front() );
   real_t factorYInner;
   uint_t indexY0;
   uint_t indexY1;
   if ( factorY >= real_c( 1.0 ) )
   {
      indexY0      = rangeY.size() - 1;
      indexY1      = indexY0;
      factorYInner = real_c( 1.0 );
   }
   else if ( factorY <= real_c( 0.0 ) )
   {
      indexY0      = 0;
      indexY1      = indexY0;
      factorYInner = real_c( 0.0 );
   }
   else
   {
      auto iY      = std::lower_bound( rangeY.begin(), rangeY.end(), y );
      indexY1      = static_cast< uint_t >( std::distance( rangeY.begin(), iY ) );
      indexY0      = indexY1 - 1;
      factorYInner = ( y - rangeY.at( indexY0 ) ) / ( rangeY.at( indexY1 ) - rangeY.at( indexY0 ) );
   }

   // get values at the corner points
   real_t f00 = values.at( indexX0 ).at( indexY0 );
   real_t f10 = values.at( indexX1 ).at( indexY0 );
   real_t f01 = values.at( indexX0 ).at( indexY1 );
   real_t f11 = values.at( indexX1 ).at( indexY1 );

   // linear interpolation in x direction
   real_t fx0 = f00 + factorXInner * ( f10 - f00 );
   real_t fx1 = f01 + factorXInner * ( f11 - f01 );

   // linear interpolation in T direction
   return fx0 + factorYInner * ( fx1 - fx0 );
}

enum NondimensionalisationType : size_t
{
   None     = 0,
   Temp     = 1,
   Space    = 2,
   Pressure = 3,
   Density  = 4,
   Cp       = 5,
   Alpha    = 6,
   KappaT   = 7,
   Gravity  = 8,
   Viscosity = 9
};

inline std::ostream& operator<<( std::ostream& os, const MantleConvection::NondimensionalisationType type )
{
   switch ( type )
   {
   case Temp:
      return os << "Temp";
   case Space:
      return os << "Space";
   case Pressure:
      return os << "Pressure";
   case Density:
      return os << "Density";
   case Cp:
      return os << "Cp";
   case Alpha:
      return os << "Alpha";
   case KappaT:
      return os << "KappaT";
   case Gravity:
      return os << "Gravity";
   case Viscosity:
      return os << "Viscosity";      
   default:
      return os << "None";
   }
}

const static std::map< std::string, MantleConvection::NondimensionalisationType > NondimensionalisationTypeMap = {
    { "Temp", MantleConvection::NondimensionalisationType::Temp },
    { "Space", MantleConvection::NondimensionalisationType::Space },
    { "None", MantleConvection::NondimensionalisationType::None },
    { "Pressure", MantleConvection::NondimensionalisationType::Pressure },
    { "Density", MantleConvection::NondimensionalisationType::Density },
    { "Cp", MantleConvection::NondimensionalisationType::Cp },
    { "Alpha", MantleConvection::NondimensionalisationType::Alpha },
    { "KappaT", MantleConvection::NondimensionalisationType::KappaT },
    { "Gravity", MantleConvection::NondimensionalisationType::Gravity },
    { "Viscosity", MantleConvection::NondimensionalisationType::Viscosity } };

struct csvData
{
   uint_t                                                                    rowsRange;
   uint_t                                                                    colsClassification;
   std::vector< std::vector< MantleConvection::NondimensionalisationType > > classification;
   std::vector< std::vector< real_t > >                                      range;
   std::vector< std::vector< real_t > >                                      values;
};

inline std::ostream& operator<<( std::ostream& os, const csvData& cd )
{
   // clang-format off
   os << "--------------csvData--------------" << "\n";
   os << "rowsRange: "           << cd.rowsRange            << "\n";
   os << "colsClassification: "  << cd.colsClassification   << "\n";
   os << "classification: "                                 << "\n";
   for (auto& vec : cd.classification)
   {
      for (auto& v : vec) {
         os << v << ", ";
      }
      if (vec.size() > 0)
      {
         os << "\n";
      }
   }
   os << "range: "                                          << "\n";
   for (auto& vec : cd.range)
   {
      for (auto& v : vec) {
         os << v << ", ";
      }
      if (vec.size() > 0)
      {
         os << "\n";
      }
   }
      os << "values: "                                      << "\n";
   for (auto& vec : cd.values)
   {
      for (auto& v : vec) {
         os << v << ", ";
      }
      if (vec.size() > 0)
      {
         os << "\n";
      }
   }
   // clang-format on

   return os;
}

void nondimensionaliseValue( NondimensionalisationParameters&             nd,
                             MantleConvection::NondimensionalisationType& type,
                             real_t&                                      value )
{
   switch ( type )
   {
   case MantleConvection::NondimensionalisationType::Temp:
      value = value / nd.deltaT_;
      break;
   case MantleConvection::NondimensionalisationType::Space:
      value = value / nd.d_;
      break;
   case MantleConvection::NondimensionalisationType::Pressure:
      value = value / nd.pRef_;
      break;
   case MantleConvection::NondimensionalisationType::Density:
      value = value / nd.rhoRef_;
      break;
   case MantleConvection::NondimensionalisationType::Cp:
      value = value / nd.C_pRef_;
      break;
   case MantleConvection::NondimensionalisationType::Alpha:
      value = value / nd.alphaRef_;
      break;
   case MantleConvection::NondimensionalisationType::KappaT:
      value = value / nd.kappa_TRef_;
      break;
   case MantleConvection::NondimensionalisationType::Gravity:
      value = value / nd.gRef_;
      break;
   case MantleConvection::NondimensionalisationType::Viscosity:
      value = value / nd.etaRef_;
      break;      
   default:
      break;
   }
}

} // namespace MantleConvection