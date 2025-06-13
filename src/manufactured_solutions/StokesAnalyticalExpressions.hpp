/*
 * Copyright (c) 2025 Marcus Mohr.
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

#include <functional>

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/math/Constants.h"

#include "hyteg/types/PointND.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::math::pi;

namespace hyteg::manufactured_solutions {

/// Enum for selecting a component expression of the Stokes problem
enum class ExpressionType
{
   VELOCITY_X,
   VELOCITY_Y,
   PRESSURE,
   RHS_X,
   RHS_Y
};

/// Map for converting component expressions to string representations
const std::map< ExpressionType, std::string > expressionType2string = { { ExpressionType::VELOCITY_X, "velocity_x" },
                                                                        { ExpressionType::VELOCITY_Y, "velocity_y" },
                                                                        { ExpressionType::PRESSURE, "pressure" },
                                                                        { ExpressionType::RHS_X, "rhs_x" },
                                                                        { ExpressionType::RHS_Y, "rhs_y" } };

/// Stationary Stokes example from the book by Volker John
///
/// This benchmark for the incompressible Stokes equation is example D.3 in the book "Finite Element Methods
/// for Incompressible Flow Problems" by Volker John. The problem is posed on the unit square and includes
/// no-slip boundary conditions for velocity.  The pressure constant is chosen such that the pressure has zero mean.
/// The analytical solution is given by
///
/// \f[
/// \begin{align*}
/// u(x,y) &= 1000 \left( \begin{array}{c}
/// x^2 (1-x)^4 y^2 (1-y) (3-5y) \\ %
/// -2x (1-x)^3 (1-3x) y^3 (1-y)^2
/// \end{array} \right) \\[2ex]
/// p(x,y) &= \pi^2 \Big[ x y^3 \cos( 2 \pi x^2 y ) - x^2 y \sin( 2 \pi x y ) \Big] + \frac{1}{8}
/// \end{align*}
/// \f]
struct Stokes_Benchmark_by_John_D3
{
   static std::function< real_t( const hyteg::Point3D& ) > getAnalyticalExpression( ExpressionType exType )
   {
      std::function< real_t( const hyteg::Point3D& ) > expression;

      switch ( exType )
      {
      case ( ExpressionType::VELOCITY_X ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            real_t y = coords[1];

            real_t tmp1 = ( real_c( 1 ) - x );
            real_t tmp2 = tmp1 * tmp1;
            real_t tmp3 = tmp2 * tmp2;

            return real_c( 1000 ) * x * x * tmp3 * y * y * ( real_c( 1 ) - y ) * ( real_c( 3 ) - real_c( 5 ) * y );
         };
         break;

      case ( ExpressionType::VELOCITY_Y ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            real_t y = coords[1];

            real_t tmp1 = ( real_c( 1 ) - x );
            real_t tmp2 = tmp1 * tmp1 * tmp1;
            real_t tmp3 = ( real_c( 1 ) - y ) * ( real_c( 1 ) - y );
            real_t y3   = y * y * y;
            return -real_c( 2000 ) * x * tmp2 * ( real_c( 1 ) - real_c( 3 ) * x ) * y3 * tmp3;
         };
         break;

      case ( ExpressionType::PRESSURE ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            real_t y = coords[1];

            real_t x2 = x * x;
            real_t y3 = y * y * y;

            return pi * pi * ( x * y3 * std::cos( real_c( 2 ) * pi * x2 * y ) - x2 * y * std::sin( real_c( 2 ) * pi * x * y ) ) +
                   real_c( 1.0 / 8.0 );
         };
         break;

      case ( ExpressionType::RHS_X ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            real_t y = coords[1];

            real_t t1  = x * x;
            real_t t4  = 0.2e1 * y * t1 * pi;
            real_t t5  = std::cos( t4 );
            real_t t6  = pi * pi;
            real_t t8  = y * y;
            real_t t9  = y * t8;
            real_t t11 = std::sin( t4 );
            real_t t12 = pi * t6;
            real_t t14 = t8 * t8;
            real_t t20 = 0.2e1 * pi * x * y;
            real_t t21 = std::cos( t20 );
            real_t t26 = std::sin( t20 );
            real_t t32 = std::pow( -0.1e1 + x, 0.2e1 );
            real_t t33 = 0.4e1 / 0.5e1 * y;
            real_t t35 = t1 * t1;
            real_t retVal =
                t9 * t6 * t5 - 0.4e1 * t14 * t1 * t12 * t11 - 0.2e1 * t8 * t1 * t12 * t21 - 0.2e1 * x * y * t6 * t26 -
                0.60000e5 *
                    ( t35 * ( t8 - t33 + 0.1e1 / 0.10e2 ) + x * t1 * ( -0.2e1 * t8 + 0.8e1 / 0.5e1 * y - 0.1e1 / 0.5e1 ) +
                      t1 * ( 0.5e1 / 0.2e1 * t14 - 0.4e1 * t9 + 0.5e1 / 0.2e1 * t8 - t33 + 0.1e1 / 0.10e2 ) +
                      x * ( -0.5e1 / 0.3e1 * t14 + 0.8e1 / 0.3e1 * t9 - t8 ) +
                      ( -0.3e1 / 0.5e1 + y ) * ( -0.1e1 + y ) * t8 / 0.6e1 ) *
                    t32;

            return retVal;
         };
         break;

      case ( ExpressionType::RHS_Y ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            real_t y = coords[1];

            real_t t1     = x * x;
            real_t t4     = 0.2e1 * y * t1 * pi;
            real_t t5     = std::cos( t4 );
            real_t t6     = pi * pi;
            real_t t8     = y * y;
            real_t t12    = std::sin( t4 );
            real_t t13    = pi * t6;
            real_t t15    = x * t1;
            real_t t16    = y * t8;
            real_t t22    = 0.2e1 * pi * x * y;
            real_t t23    = std::cos( t22 );
            real_t t28    = std::sin( t22 );
            real_t t31    = -0.1e1 + x;
            real_t t34    = t8 * t8;
            real_t t40    = t1 * t1;
            real_t t46    = t31 * t31;
            real_t t48    = ( -0.1e1 / 0.3e1 + x ) * t46;
            real_t retVal = 0.3e1 * t8 * x * t6 * t5 - 0.2e1 * t16 * t15 * t13 * t12 - 0.2e1 * y * t15 * t13 * t23 -
                            t1 * t6 * t28 +
                            0.120000e6 *
                                ( t34 * ( t1 - x + 0.1e1 / 0.5e1 ) + t16 * ( -0.2e1 * t1 + 0.2e1 * x - 0.2e1 / 0.5e1 ) +
                                  t8 * ( t40 - 0.7e1 / 0.3e1 * t15 + 0.8e1 / 0.3e1 * t1 - 0.4e1 / 0.3e1 * x + 0.1e1 / 0.5e1 ) -
                                  0.6e1 / 0.5e1 * x * y * t48 + 0.3e1 / 0.10e2 * x * t48 ) *
                                y * t31;

            return retVal;
         };
         break;
      }

      return expression;
   }
};

/// A Stokes example from the book by Jean Donea and Antonio Huerta
///
/// This is a manufactured example of a simple stationary 2D convection problem. It is taken from the book
/// "Finite Element Methods for Flow Problems", John Wiley & Sons (2003) by Jean Donea and Antonio Huerta.
///
/// The flow is incompressible and isoviscous, with both viscosity and density taken to be one. Hence the
/// forcing term is simply given by the gravity field. The latter is artificially chosen to get a rotating
/// flow field. The problem is posed on the unit square and assumes no-slip boundary conditions for velocity. The pressure constant is chosen such that the pressure has zero mean.
///
/// The forcing term and analytical solution are given by
///
/// \f[
/// \begin{align*}
/// u(x,y) &= \left( \begin{array}{ll} x^2(1- x)^2 (2y - 6y^2 + 4y^3) \\ %
/// -y^2 (1 - y)^2 (2x - 6x^2 + 4x^3) \end{array} \right) \\[2ex]
/// p(x,y) &=  x(1 -x) - \frac{1}{6}\\[2ex]
/// g(x,y) &= \left( \begin{array}{c}
/// (12 - 24y) x^4 + (-24 + 48y) x^3 + (-48y + 72y^2 - 48 y^3 + 12) x^2 + (-2 + 24y -72y^2+48y^3)x + 1-4y + 12y^2-8y^3 \\
/// (8 - 48y + 48 y^2) x^3 + (-12 + 72y - 72y^2) x^2 + (4 - 24y + 48y^2 - 48y^3 + 24y^4) x - 12y^2 + 24y^3 - 12y^4
/// \end{array} \right)
/// \end{align*}
/// \f]
///
/// \note This is one of the benchmarks used in "On the choice of finite element for applications in geodynamics"
/// by Cedric Thieulot1 and Wolfgang Bangerth, <a href="https://doi.org/10.5194/se-13-229-2022">10.5194/se-13-229-2022</a>
struct Stokes_Benchmark_by_Donea_Huerta
{
   static std::function< real_t( const hyteg::Point3D& ) > getAnalyticalExpression( ExpressionType exType )
   {
      std::function< real_t( const hyteg::Point3D& ) > expression;
      switch ( exType )
      {
      case ( ExpressionType::VELOCITY_X ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            real_t y = coords[1];

            real_t tmp1 = x * x;
            real_t tmp2 = real_c( 1 ) - x;
            real_t tmp3 = y * y;
            real_t tmp4 = tmp3 * y;

            return tmp1 * tmp2 * tmp2 * ( real_c( 2 ) * y - real_c( 6 ) * tmp3 + real_c( 4 ) * tmp4 );
         };
         break;

      case ( ExpressionType::VELOCITY_Y ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            real_t y = coords[1];

            real_t tmp1 = y * y;
            real_t tmp2 = real_c( 1 ) - y;
            real_t tmp3 = x * x;
            real_t tmp4 = tmp3 * x;

            return -tmp1 * tmp2 * tmp2 * ( real_c( 2 ) * x - real_c( 6 ) * tmp3 + real_c( 4 ) * tmp4 );
         };
         break;

      case ( ExpressionType::PRESSURE ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            return x * ( real_c( 1 ) - x ) - real_c( 1.0 / 6.0 );
         };
         break;

      case ( ExpressionType::RHS_X ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            real_t y = coords[1];

            real_t t2  = 0.1e1 - x;
            real_t t3  = t2 * t2;
            real_t t4  = y * y;
            real_t t9  = 0.4e1 * t4 * y - 0.6e1 * t4 + 0.2e1 * y;
            real_t t15 = x * x;

            return 0.1e1 - 0.2e1 * x - 0.2e1 * t9 * t3 + 0.8e1 * t9 * t2 * x - 0.2e1 * t9 * t15 -
                   ( 0.24e2 * y - 0.12e2 ) * t3 * t15;
         };
         break;

      case ( ExpressionType::RHS_Y ):
         expression = []( const hyteg::Point3D& coords ) {
            real_t x = coords[0];
            real_t y = coords[1];

            real_t t1  = y * y;
            real_t t2  = 0.1e1 - y;
            real_t t3  = t2 * t2;
            real_t t8  = x * x;
            real_t t13 = 0.4e1 * t8 * x - 0.6e1 * t8 + 0.2e1 * x;

            return ( 0.24e2 * x - 0.12e2 ) * t3 * t1 + 0.2e1 * t13 * t3 - 0.8e1 * t13 * t2 * y + 0.2e1 * t13 * t1;
         };
         break;
      }

      return expression;
   }
};

/// Accessor class for obtaining analytical expressions for different 2D Stokes benchmarks
struct Stokes2D
{
   enum class BenchmarkType
   {
      JOHN_D3,
      DONEA_HUERTA
   };

   /// Get an analytical function description of a benchmark component
   ///
   /// Returns an expression in the form of an std::function object that allows to evaluate the analytic solution or the
   /// right-hand side of the requested 2D Stokes benchmark
   static std::function< real_t( const hyteg::Point3D& ) > get( Stokes2D::BenchmarkType benchmarkType, ExpressionType exType )
   {
      switch ( benchmarkType )
      {
      case ( Stokes2D::BenchmarkType::JOHN_D3 ):
         return Stokes_Benchmark_by_John_D3::getAnalyticalExpression( exType );

      case ( Stokes2D::BenchmarkType::DONEA_HUERTA ):
         return Stokes_Benchmark_by_Donea_Huerta::getAnalyticalExpression( exType );

      default:
         WALBERLA_ABORT( "Ooops. Unhandled benchmarkType!" );
      }
   }
};

} // namespace hyteg::manufactured_solutions
