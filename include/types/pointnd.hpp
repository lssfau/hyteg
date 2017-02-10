#ifndef POINT3D_HPP
#define POINT3D_HPP

#include <fmt/ostream.h>
#include <cmath>
#include <array>

namespace hhg
{

template<typename T, size_t N>
class PointND
{

public:
  PointND()
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] = (T) 0;
    }
  }

  PointND(T _x[N])
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] = _x[i];
    }
  }

  PointND(std::array<T, N> list)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] = list[i];
    }
  }

  PointND(const PointND& b)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] = b.x[i];
    }
  }

  T dot(const PointND& b) const
  {
    T tmp = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
      tmp += x[i] * b.x[i];
    }
    return tmp;
  }

  T normSq() const
  {
    return dot(*this);
  }

  T norm() const
  {
    return std::sqrt(normSq());
  }

  PointND& operator+=(const PointND& rhs)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] += rhs.x[i];
    }
    return *this;
  }

  PointND& operator-=(const PointND& rhs)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] -= rhs.x[i];
    }
    return *this;
  }

  PointND& operator*=(T scalar)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] *= scalar;
    }
    return *this;
  }

  PointND& operator/=(T scalar)
  {
    for (size_t i = 0; i < N; ++i)
    {
      x[i] /= scalar;
    }
    return *this;
  }

  T& operator[](const int index)
  {
    return x[index];
  }

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

typedef PointND<double, 2> Point2D;
typedef PointND<double, 3> Point3D;

}

#endif /* POINT3D_HPP */