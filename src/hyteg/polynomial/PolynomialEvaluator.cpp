/*
* Copyright (c) 2017-2024 Dominik Thoennes
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

#include "PolynomialEvaluator.hpp"

namespace hyteg {
void Polynomial2DEvaluator::setY( real_t y )
{
   // std::cout << "Evaluator2D(q=" << degree_ << ")::setY of p2 with q = "<< poly2_->getDegree() << "\n";

   for ( uint_t degree = 0; degree <= degree_; ++degree )
   {
      poly1_.setCoefficient( degree, 0.0 );
   }

   uint_t start = 0;
   real_t y_;

   for ( uint_t coeff = 0; coeff <= degree_; ++coeff )
   {
      uint_t idx = start;
      y_         = real_t( 1.0 );

      for ( uint_t degree = 0; degree <= degree_ - coeff; ++degree )
      {
         poly1_.addToCoefficient( coeff, poly2_->getCoefficient( idx ) * y_ );

         idx += coeff + degree + 2;
         y_ *= y;
      }

      start += coeff + 1;
   }
}

real_t Polynomial2DEvaluator::evalX( real_t x ) const
{
   // return poly1_.eval(x);
   real_t px = poly1_.getCoefficient( 0 );
   real_t xk = 1.0;
   for ( uint_t k = 1; k <= degree_; ++k )
   {
      xk *= x;
      px += poly1_.getCoefficient( k ) * xk;
   }
   return px;
}

real_t Polynomial2DEvaluator::setStartX( real_t x, real_t h )
{
   switch ( degree_ )
   {
   case 0:
      return polynomialevaluator::setStartX< 0 >( x, h, poly1_, deltas_ );
   case 1:
      return polynomialevaluator::setStartX< 1 >( x, h, poly1_, deltas_ );
   case 2:
      return polynomialevaluator::setStartX< 2 >( x, h, poly1_, deltas_ );
   case 3:
      return polynomialevaluator::setStartX< 3 >( x, h, poly1_, deltas_ );
   case 4:
      return polynomialevaluator::setStartX< 4 >( x, h, poly1_, deltas_ );
   case 5:
      return polynomialevaluator::setStartX< 5 >( x, h, poly1_, deltas_ );
   case 6:
      return polynomialevaluator::setStartX< 6 >( x, h, poly1_, deltas_ );
   case 7:
      return polynomialevaluator::setStartX< 7 >( x, h, poly1_, deltas_ );
   case 8:
      return polynomialevaluator::setStartX< 8 >( x, h, poly1_, deltas_ );
   case 9:
      return polynomialevaluator::setStartX< 9 >( x, h, poly1_, deltas_ );
   case 10:
      return polynomialevaluator::setStartX< 10 >( x, h, poly1_, deltas_ );
   case 11:
      return polynomialevaluator::setStartX< 11 >( x, h, poly1_, deltas_ );
   case 12:
      return polynomialevaluator::setStartX< 12 >( x, h, poly1_, deltas_ );
   default:
      return 0;
   }
}

real_t Polynomial2DEvaluator::incrementEval()
{
   switch ( degree_ )
   {
   case 0:
      return polynomialevaluator::incrementEval< 0 >( deltas_ );
   case 1:
      return polynomialevaluator::incrementEval< 1 >( deltas_ );
   case 2:
      return polynomialevaluator::incrementEval< 2 >( deltas_ );
   case 3:
      return polynomialevaluator::incrementEval< 3 >( deltas_ );
   case 4:
      return polynomialevaluator::incrementEval< 4 >( deltas_ );
   case 5:
      return polynomialevaluator::incrementEval< 5 >( deltas_ );
   case 6:
      return polynomialevaluator::incrementEval< 6 >( deltas_ );
   case 7:
      return polynomialevaluator::incrementEval< 7 >( deltas_ );
   case 8:
      return polynomialevaluator::incrementEval< 8 >( deltas_ );
   case 9:
      return polynomialevaluator::incrementEval< 9 >( deltas_ );
   case 10:
      return polynomialevaluator::incrementEval< 10 >( deltas_ );
   case 11:
      return polynomialevaluator::incrementEval< 11 >( deltas_ );
   case 12:
      return polynomialevaluator::incrementEval< 12 >( deltas_ );
   default:
      return 0;
   }
}

void Polynomial3DEvaluator::setZ( real_t z )
{
   poly_z_.setZero();

   // idx of coefficient c_ijk of 3d polynoial
   int idx_ijk = 0;

   // p(x,y,z) = sum_{ |i+j+k| = 0,...,degree_ } c_ijk * x^i*y^j*z^k
   for ( int ijk = 0; ijk <= degree_; ++ijk )
   {
      for ( int i = ijk; i >= 0; --i )
      {
         // idx of coefficient c_ij of 2d polynoial
         int idx_ij = ijk * ( ijk + 1 ) / 2 + ijk - i;
         // z^k
         real_t z_k = real_t( 1.0 );

         /* compute coefficients of 2d polynomial:
          p|_{z}(x,y) = sum{ |i+j| = 0,...,degree_ } ( sum_{ k = 0,...,degree_ - |i+j| } c_ijk * z^k ) * x^i*y^j
        */
         for ( int k = 0; k <= ijk - i; ++k ) // note: j = ijk - i - k
         {
            // c_ij += c_ijk * z^k
            poly_z_.addToCoefficient( walberla::uint_c( idx_ij ), poly_->getCoefficient( walberla::uint_c( idx_ijk ) ) * z_k );

            idx_ij -= ( ijk - k + 1 );
            ++idx_ijk;
            z_k *= z;
         }
      }
   }
}

void Polynomial3DEvaluator::setY( real_t y )
{
   poly_yz_.setZero();

   uint_t start = 0;
   real_t y_k;

   for ( uint_t coeff = 0; coeff <= degree_; ++coeff )
   {
      uint_t idx = start;
      y_k        = real_t( 1.0 );

      for ( uint_t degree = 0; degree <= degree_ - coeff; ++degree )
      {
         poly_yz_.addToCoefficient( coeff, poly_z_.getCoefficient( idx ) * y_k );

         idx += coeff + degree + 2;
         y_k *= y;
      }

      start += coeff + 1;
   }
}

real_t Polynomial3DEvaluator::evalX( real_t x ) const
{
   // return poly_yz_.eval(x);
   real_t px = poly_yz_.getCoefficient( 0 );
   real_t xk = 1.0;
   for ( uint_t k = 1; k <= degree_; ++k )
   {
      xk *= x;
      px += poly_yz_.getCoefficient( k ) * xk;
   }
   return px;
}

real_t Polynomial3DEvaluator::setStartX( real_t x, real_t h )
{
   switch ( degree_ )
   {
   case 0:
      return polynomialevaluator::setStartX< 0 >( x, h, poly_yz_, deltas_ );
   case 1:
      return polynomialevaluator::setStartX< 1 >( x, h, poly_yz_, deltas_ );
   case 2:
      return polynomialevaluator::setStartX< 2 >( x, h, poly_yz_, deltas_ );
   case 3:
      return polynomialevaluator::setStartX< 3 >( x, h, poly_yz_, deltas_ );
   case 4:
      return polynomialevaluator::setStartX< 4 >( x, h, poly_yz_, deltas_ );
   case 5:
      return polynomialevaluator::setStartX< 5 >( x, h, poly_yz_, deltas_ );
   case 6:
      return polynomialevaluator::setStartX< 6 >( x, h, poly_yz_, deltas_ );
   case 7:
      return polynomialevaluator::setStartX< 7 >( x, h, poly_yz_, deltas_ );
   case 8:
      return polynomialevaluator::setStartX< 8 >( x, h, poly_yz_, deltas_ );
   case 9:
      return polynomialevaluator::setStartX< 9 >( x, h, poly_yz_, deltas_ );
   case 10:
      return polynomialevaluator::setStartX< 10 >( x, h, poly_yz_, deltas_ );
   case 11:
      return polynomialevaluator::setStartX< 11 >( x, h, poly_yz_, deltas_ );
   case 12:
      return polynomialevaluator::setStartX< 12 >( x, h, poly_yz_, deltas_ );
   default:
      return 0;
   }
}

real_t Polynomial3DEvaluator::incrementEval()
{
   switch ( degree_ )
   {
   case 0:
      return polynomialevaluator::incrementEval< 0 >( deltas_ );
   case 1:
      return polynomialevaluator::incrementEval< 1 >( deltas_ );
   case 2:
      return polynomialevaluator::incrementEval< 2 >( deltas_ );
   case 3:
      return polynomialevaluator::incrementEval< 3 >( deltas_ );
   case 4:
      return polynomialevaluator::incrementEval< 4 >( deltas_ );
   case 5:
      return polynomialevaluator::incrementEval< 5 >( deltas_ );
   case 6:
      return polynomialevaluator::incrementEval< 6 >( deltas_ );
   case 7:
      return polynomialevaluator::incrementEval< 7 >( deltas_ );
   case 8:
      return polynomialevaluator::incrementEval< 8 >( deltas_ );
   case 9:
      return polynomialevaluator::incrementEval< 9 >( deltas_ );
   case 10:
      return polynomialevaluator::incrementEval< 10 >( deltas_ );
   case 11:
      return polynomialevaluator::incrementEval< 11 >( deltas_ );
   case 12:
      return polynomialevaluator::incrementEval< 12 >( deltas_ );
   default:
      return 0;
   }
}

} // namespace hyteg
