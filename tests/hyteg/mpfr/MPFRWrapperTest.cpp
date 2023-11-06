/*
* Copyright (c) 2023-2024 Michael Zikeli.
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

#include <inttypes.h>
#include <string_view>
#include <sstream>
#include <stdio.h>

#ifndef MPFR_USE_INTMAX_T
#define MPFR_USE_INTMAX_T
#endif //MPFR_USE_INTMAX_T
#ifndef MPFR_USE_NO_MACRO
#define MPFR_USE_NO_MACRO
#endif //MPFR_USE_NO_MACRO

#include <mpfr.h>

#include "core/Environment.h"
#include "core/debug/Debug.h"

#include "mpfr_real_v0.0.9-alpha.hpp"
#include "mpreal/mpreal.h"

namespace hyteg::mpfrTest {
constexpr unsigned int mpfrAccuracy = 200;
constexpr unsigned int ouputMantissa = 61;
constexpr std::string_view exactSolutionString = "2.7182818284590452353602874713526624977572470936999595749669704e0";
// NOTE: To get the exact solution, the exact example from the MPFR website was used.
// After the correctness of the 'checkPureMPFR()' function was verified, the result was copied is now this string_view.

// This example comes from the MPFR website.
// Only the rounding is altered to "round to nearest" (MPFR_RNDN) since it's the standard for both wrapper.
// Compare https://www.mpfr.org/sample.html
void checkPureMPFR()
{
   unsigned int i;
   mpfr_t       s, t, u;

   mpfr_init2( t, mpfrAccuracy );
   mpfr_set_d( t, 1.0, MPFR_RNDN );
   mpfr_init2( s, mpfrAccuracy );
   mpfr_set_d( s, 1.0, MPFR_RNDN );
   mpfr_init2( u, mpfrAccuracy );
   for ( i = 1; i <= 100; i++ )
   {
      mpfr_mul_ui( t, t, i, MPFR_RNDN );
      mpfr_set_d( u, 1.0, MPFR_RNDN );
      mpfr_div( u, u, t, MPFR_RNDN );
      mpfr_add( s, s, u, MPFR_RNDN );
   }

   printf( "Sum is  " );
   mpfr_out_str( stdout, 10, 0, s, MPFR_RNDN );
   putchar( '\n' );

   mpfr_t exact;           // The exact solution is from the MPFR tutorial from the website
   mpfr_t smaller, larger; // smaller and larger are chosen arbitrarily
   mpfr_init2( exact, mpfrAccuracy );
   mpfr_set_str( exact, exactSolutionString.data(), 0, MPFR_RNDN );
   mpfr_init2( smaller, mpfrAccuracy );
   mpfr_set_str( smaller, "-2.7182818284590452353602874713526624977572470936999595749669704e0", 0, MPFR_RNDN );
   mpfr_init2( larger, mpfrAccuracy );
   mpfr_set_str( larger, "2.7182818284590452353602874713526624977572470936999595749669724e0", 0, MPFR_RNDN );

   printf( "Exa is  " );
   mpfr_out_str( stdout, 10, 0, exact, MPFR_RNDN );
   putchar( '\n' );
   printf( "Low is " );
   mpfr_out_str( stdout, 10, 0, smaller, MPFR_RNDN );
   putchar( '\n' );
   printf( "Hig is  " );
   mpfr_out_str( stdout, 10, 0, larger, MPFR_RNDN );
   putchar( '\n' );

   WALBERLA_CHECK_EQUAL( mpfr_cmp( exact, s ), 0 );
   WALBERLA_CHECK_LESS( mpfr_cmp( smaller, s ), 0 );
   WALBERLA_CHECK_GREATER( mpfr_cmp( larger, s ), 0 );

   mpfr_clear( s );
   mpfr_clear( t );
   mpfr_clear( u );
   mpfr_clear( exact );
   mpfr_free_cache();

   return;
} // checkPureMPFR()

void compareWithCompileTimeWrapper()
{
   unsigned int               i;
   mpfr::real< mpfrAccuracy > s = "1.0", t = "1.0", u;

   for ( i = 1; i <= 100; i++ )
   {
      mpfr::real< mpfrAccuracy > tmp = i;
      t *= tmp;
      u = "1.0";
      u /= t;
      s += u;
   }

   std::cout << "compile-time wrapper : " << s << std::endl;
   std::cout << "Exact                : " << exactSolutionString.data() << std::endl;

   mpfr_t exact;           // The exact solution is from the MPFR tutorial from the website
   mpfr_init2( exact, mpfrAccuracy );
   mpfr_set_str( exact, exactSolutionString.data(), 0, MPFR_RNDN );
   WALBERLA_CHECK_EQUAL( mpfr::cmpabs( exact, s ), 0 );

   return;
} // compareWithCompileTimeWrapper()

void compareWithRuntimeWrapper()
{
   unsigned int i;
   mpfr::mpreal s("1.0", mpfrAccuracy), t("1.0", mpfrAccuracy), u(0, mpfrAccuracy);

   for ( i = 1; i <= 100; i++ )
   {
      t *= i;
      u = "1.0";
      u /= t;
      s += u;
   }
   std::cout << "Runtime wrapper : " << s << std::endl;
   std::cout << "Exact           : " << exactSolutionString.data() << std::endl;

   mpfr_t exact;           // The exact solution is from the MPFR tutorial from the website
   mpfr_init2( exact, mpfrAccuracy );
   mpfr_set_str( exact, exactSolutionString.data(), 0, MPFR_RNDN );
   WALBERLA_CHECK_EQUAL( mpfr::cmpabs( exact, s ), 0 );

   return;
} // compareWithRuntimeWrapper()

} // namespace hyteg::mpfrTest

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   std::cout << std::setbase(10) << std::scientific << std::setprecision(hyteg::mpfrTest::ouputMantissa);
   hyteg::mpfrTest::checkPureMPFR();
   hyteg::mpfrTest::compareWithCompileTimeWrapper();
   hyteg::mpfrTest::compareWithRuntimeWrapper();

   return EXIT_SUCCESS;
}
