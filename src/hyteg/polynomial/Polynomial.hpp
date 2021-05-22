/*
 * Copyright (c) 2021 Daniel Drzisga, Dominik Thoennes, Benjamin Mann
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

#include "MonomialBasis1D.hpp"
#include "MonomialBasis2D.hpp"
#include "MonomialBasis3D.hpp"

#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

namespace hyteg {

using walberla::mpi::RecvBuffer;
using walberla::mpi::SendBuffer;


template <int Dim, class Point, class Basis>
class Polynomial
{
 public:
   static constexpr uint_t getNumCoefficientsForDegree(uint_t degree)
   {
      return math::binomialCoefficient(Dim + degree - 1, degree);
   }

   static constexpr uint_t getNumCoefficients(uint_t degree)
   {
      return math::binomialCoefficient(Dim + degree, degree);
   }

   Polynomial(uint_t degree = 0)
      : degree_(degree),
        numCoefficients_(getNumCoefficients(degree)),
        coeffs_(numCoefficients_)
   {
   }

   Polynomial(Polynomial&& p) : Polynomial(p.getDegree())
   {
      coeffs_ = p.coeffs_;
   }

   Polynomial(const Polynomial& p) : Polynomial(p.getDegree())
   {
      coeffs_ = p.coeffs_;
   }

   inline Polynomial& operator=(Polynomial&& p)
   {
      WALBERLA_ASSERT(degree_ == p.getDegree(), "Polynomial degrees don't match!");
      coeffs_ = p.coeffs_;
      return *this;
   }

   inline Polynomial& operator=(const Polynomial& p)
   {
      WALBERLA_ASSERT(degree_ == p.getDegree(), "Polynomial degrees don't match!");
      coeffs_ = p.coeffs_;
      return *this;
   }

   inline uint_t getDegree() const
   {
      return degree_;
   }

   inline real_t eval(const Point& x) const
   {
      real_t eval = coeffs_[0] * Basis::eval(0, x);

      for (uint_t c = 1; c < numCoefficients_; ++c)
      {
         eval = std::fma(coeffs_[c], Basis::eval(c, x), eval);
      }

      return eval;
   }

   inline void setCoefficient(uint_t idx, real_t value)
   {
      WALBERLA_ASSERT(idx < numCoefficients_);
      coeffs_[idx] = value;
   }

   inline real_t getCoefficient(uint_t idx) const
   {
      WALBERLA_ASSERT(idx < numCoefficients_);
      return coeffs_[idx];
   }

   inline void addToCoefficient(uint_t idx, real_t value)
   {
      WALBERLA_ASSERT(idx < numCoefficients_);
      coeffs_[idx] += value;
   }

   inline void scale(real_t scalar)
   {
      for (uint_t i = 0; i < numCoefficients_; ++i)
      {
         coeffs_[i] *= scalar;
      }
   }

   inline void scaleAdd(real_t scalar, const Polynomial& rhs)
   {
      WALBERLA_ASSERT(degree_ == rhs.degree_);

      for (uint_t i = 0; i < numCoefficients_; ++i)
      {
         coeffs_[i] += scalar * rhs.coeffs_[i];
      }
   }

   inline void setZero()
   {
      std::memset(coeffs_.data(), 0, numCoefficients_ * sizeof(real_t));
   }

   real_t lInfinityError(const Polynomial& rhs)
   {
      real_t error = std::numeric_limits<real_t>::min();

      uint_t i = 0;

      for (; i < std::min(numCoefficients_, rhs.numCoefficients_); ++i)
      {
         error = std::max(error, std::abs(coeffs_[i] - rhs.coeffs_[i]));
      }

      for (; i < numCoefficients_; ++i)
      {
         error = std::max(error, std::abs(coeffs_[i]));
      }

      for (; i < rhs.numCoefficients_; ++i)
      {
         error = std::max(error, std::abs(rhs.coeffs_[i]));
      }

      return error;
   }

   /// Serializes the allocated data to a send buffer
   inline void serialize(SendBuffer& sendBuffer) const
   {
      sendBuffer << degree_;
      sendBuffer << numCoefficients_;
      sendBuffer << coeffs_;
   }

   /// Deserializes data from a recv buffer (clears all already allocated data and replaces it with the recv buffer's content)
   inline void deserialize(RecvBuffer& recvBuffer)
   {
      coeffs_.clear();

      recvBuffer >> degree_;
      recvBuffer >> numCoefficients_;
      recvBuffer >> coeffs_;
   }


 private:
   uint_t degree_;
   uint_t numCoefficients_;
   std::vector<real_t> coeffs_;

};

template <int Dim, class Point, class Basis>
inline std::ostream& operator<<(std::ostream& os, const Polynomial<Dim, Point, Basis>& poly)
{
   os << "[" << poly.getCoefficients(0);

   for (size_t i = 1; i < poly.getNumCoefficients(poly.getDegree()); ++i)
   {
      os << ", " << poly.getCoefficient(i);
   }

   os << "]";

   return os;
}

template<class Basis>
using Polynomial1D = Polynomial<1, real_t,  Basis>;
template<class Basis>
using Polynomial2D = Polynomial<2, Point2D, Basis>;
template<class Basis>
using Polynomial3D = Polynomial<3, Point3D, Basis>;

using GeneralPolynomial1D = Polynomial1D<MonomialBasis1D>;
using GeneralPolynomial2D = Polynomial2D<MonomialBasis2D>;
using GeneralPolynomial3D = Polynomial3D<MonomialBasis3D>;

} // namespace hyteg

namespace walberla {
namespace mpi {

template < typename T, // Element type of SendBuffer
           typename G, // Growth policy of SendBuffer
           int Dim, class Point, class Basis>
inline mpi::GenericSendBuffer< T, G >& operator<<(
   mpi::GenericSendBuffer< T, G >& buffer, const hyteg::Polynomial<Dim, Point, Basis>& p)
{
   p.serialize(buffer);
   return buffer;
}

template < typename T, // Element type  of RecvBuffer
           int Dim, class Point, class Basis>
inline mpi::GenericRecvBuffer< T >& operator>>(
   mpi::GenericRecvBuffer< T >& buffer, hyteg::Polynomial<Dim, Point, Basis>& p)
{
   p.deserialize(buffer);
   return buffer;
}

template <int Dim, class Point, class Basis>
struct BufferSizeTrait< hyteg::Polynomial<Dim, Point, Basis>>
{
   static const bool constantSize = false;
};

} // namespace mpi
} // namespace walberla
