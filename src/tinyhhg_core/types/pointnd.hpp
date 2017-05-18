#ifndef POINT3D_HPP
#define POINT3D_HPP

#include "core/DataTypes.h"

#include <fmt/ostream.h>
#include <cmath>
#include <array>

namespace hhg
{

using walberla::real_t;

/// \brief  N-dimensional vector
/// \author Daniel Drzisga (drzisga@ma.tum.de)
/// \date   March, 2017
///
/// The PointND class represents an N-dimensional vector with support for algebraic operations
/// \tparam T Vector data type
/// \tparam N Dimension of vector
template<typename T, size_t N>
class PointND
{

public:
  /// Constructor setting all components to zero
  PointND()
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] = (T) 0;
    }
  }

  /// Constructs the vector using values from n-dimensional array \p _x
  /// \param _x Pointer to N-dimensional array
  PointND(T _x[N])
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] = _x[i];
    }
  }

  /// Constructs the vector using values from n-dimensional array \p list, required for list initializer
  /// \param list N-dimensional array
  PointND(std::array<T, N> list)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] = list[i];
    }
  }

  /// Copy constructor
  /// \param b Reference to another instance of PointND
  PointND(const PointND& b)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] = b.x[i];
    }
  }

  /// Computes the dot product between two PointND vectors
  /// \param b Right hand side of dot operator
  /// \returns Dot product between \p this and \p b
  T dot(const PointND& b) const
  {
    T tmp = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
      tmp += x[i] * b.x[i];
    }
    return tmp;
  }

  /// Computes the squared Euclidean norm of \p this
  /// \returns Squared Euclidean norm of \p this
  T normSq() const
  {
    return dot(*this);
  }

  /// Computes the Euclidean norm of \p this
  /// \returns Euclidean norm of \p this
  T norm() const
  {
    return std::sqrt(normSq());
  }

  /// Add another PointND component wise to \p this
  /// \param rhs Reference to PointND that will be added to \p this
  /// \returns Reference to \p this after addition
  PointND& operator+=(const PointND& rhs)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] += rhs.x[i];
    }
    return *this;
  }

  /// Subtract another PointND component wise from \p this
  /// \param rhs Reference to PointND that will be subtracted from \p this
  /// \returns Reference to \p this after subtraction
  PointND& operator-=(const PointND& rhs)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] -= rhs.x[i];
    }
    return *this;
  }

  /// Multiply \p this with scalar value
  /// \param scalar Scalar value that \p this gets multiplied with
  /// \returns Reference to \p this after multiplication with \p scalar
  PointND& operator*=(T scalar)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] *= scalar;
    }
    return *this;
  }

  /// Divide \p this with scalar value
  /// \param scalar Scalar value that \p gets divided by
  /// \returns Reference to \p this after division by \p scalar
  PointND& operator/=(T scalar)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] /= scalar;
    }
    return *this;
  }

  /// Reference to component of vector at position \p index
  /// \param index The index of the component to access
  /// \returns Reference to component at position \p index
  T& operator[](const int index)
  {
    return x[index];
  }

  /// Value of component of vector at position \p index
  /// \param index The index of the component to access
  /// \returns Value of component at position \p index
  T operator[](const int index) const
  {
    return x[index];
  }

  T x[N];
};

template<typename T, size_t N>
inline PointND<T, N> operator+(PointND<T, N> lhs, const PointND<T, N>& rhs)
{
  return lhs += rhs;
}

template<typename T, size_t N>
inline PointND<T, N> operator-(PointND<T, N> lhs, const PointND<T, N>& rhs)
{
  return lhs -= rhs;
}

template<typename T, size_t N>
inline PointND<T, N> operator*(double scalar, PointND<T, N> rhs)
{
  rhs *= scalar;
  return rhs;
}

template<typename T, size_t N>
inline PointND<T, N> operator*(PointND<T, N> lhs, double scalar)
{
  lhs *= scalar;
  return lhs;
}

template<typename T, size_t N>
inline PointND<T, N> operator/(PointND<T, N> lhs, double scalar)
{
  lhs /= scalar;
  return lhs;
}

template<typename T, size_t N>
inline std::ostream& operator<<(std::ostream &os, const PointND<T, N> &pointnd)
{
  os << "[";

  for (size_t i = 0; i < N; ++i)
  {
    os << pointnd[i];
    if (i != N-1)
    {
      os << ", ";
    }
  }

  os << "]";

  return os;
}

typedef PointND<real_t, 2> Point2D;
typedef PointND<real_t, 3> Point3D;

}

#endif /* POINT3D_HPP */