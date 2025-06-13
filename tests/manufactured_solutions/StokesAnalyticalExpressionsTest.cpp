/*
 * Copyright (c) 2025 Marcus Mohr
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
#include "manufactured_solutions/StokesAnalyticalExpressions.hpp"

#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg::manufactured_solutions;

void test_John_D3()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Checking Benchmark Expressions for 'John_D3' for regressions:" );

   const std::map< ExpressionType, real_t > control = { { ExpressionType::VELOCITY_X, real_c( 0.2232790312 ) },
                                                        { ExpressionType::VELOCITY_Y, real_c( -3.438497082 ) },
                                                        { ExpressionType::PRESSURE, real_c( 0.2707198531 ) },
                                                        { ExpressionType::RHS_X, real_c( 7.585559187 ) },
                                                        { ExpressionType::RHS_Y, real_c( -463.9720348 ) } };

   for ( const auto& iter : expressionType2string )
   {
      auto   expression = Stokes2D::get( Stokes2D::BenchmarkType::JOHN_D3, iter.first );
      real_t value      = expression( { real_c( 0.10 ), real_c( 0.55 ), real_c( 0.00 ) } );
      WALBERLA_LOG_INFO_ON_ROOT( "-> " << iter.second << "(0.10,0.55) evaluates to " << std::scientific << std::setprecision( 12 )
                                       << value );
      WALBERLA_CHECK_FLOAT_EQUAL( control.at( iter.first ), value );
   }
}

void test_Donea_Huerta()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Checking Benchmark Expressions for 'Donea_Huerta' for regressions:" );

   const std::map< ExpressionType, real_t > control = { { ExpressionType::VELOCITY_X, real_c( 6.591796875000e-03 ) },
                                                        { ExpressionType::VELOCITY_Y, real_c( -6.591796875000e-03 ) },
                                                        { ExpressionType::PRESSURE, real_c( 2.083333333333e-02 ) },
                                                        { ExpressionType::RHS_X, real_c( 7.578125000000e-01 ) },
                                                        { ExpressionType::RHS_Y, real_c( -2.578125000000e-01 ) } };

   for ( const auto& iter : expressionType2string )
   {
      auto   expression = Stokes2D::get( Stokes2D::BenchmarkType::DONEA_HUERTA, iter.first );
      real_t value      = expression( { real_c( 0.25 ), real_c( 0.25 ), real_c( 0.00 ) } );
      WALBERLA_LOG_INFO_ON_ROOT( "-> " << iter.second << "(0.25,0.25) evaluates to " << std::scientific << std::setprecision( 12 )
                                       << value );
      WALBERLA_CHECK_FLOAT_EQUAL( control.at( iter.first ), value );
   }
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   test_Donea_Huerta();
   test_John_D3();
}
