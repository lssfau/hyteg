#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/communication/BufferedCommunication.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"


namespace hhg {



template< typename ValueType >
class P2Function : public Function< P2Function< ValueType > >
{
public:

  P2Function( const std::string& name, const std::shared_ptr< PrimitiveStorage > & storage, uint_t minLevel, uint_t maxLevel ) :
      Function< P2Function< ValueType > >( name, storage, minLevel, maxLevel ),
      vertexDoFFunction(name + "_v", storage, minLevel, maxLevel),
      edgeDoFFunction(name + "_e", storage, minLevel, maxLevel)
  {
  }

  P1Function<ValueType>& getVertexDoFFunction() {
    return vertexDoFFunction;
  }

  const P1Function<ValueType>& getVertexDoFFunction() const {
    return vertexDoFFunction;
  }

  EdgeDoFFunction<ValueType>& getEdgeDoFFunction() {
    return edgeDoFFunction;
  }

  const EdgeDoFFunction<ValueType>& getEdgeDoFFunction() const {
    return edgeDoFFunction;
  }

private:

  P1Function<ValueType> vertexDoFFunction;
  EdgeDoFFunction<ValueType> edgeDoFFunction;

  /// Interpolates a given expression to a P2Function
  inline void interpolate_impl(std::function<ValueType(const Point3D&)>& expr, uint_t level, DoFType flag = All);

  inline void assign_impl(const std::vector<ValueType> scalars, const std::vector<P2Function< ValueType >*> functions, uint_t level, DoFType flag = All);

  inline void add_impl(const std::vector<ValueType> scalars, const std::vector<P2Function< ValueType >*> functions, uint_t level, DoFType flag = All);

  inline real_t dot_impl(P2Function< ValueType >& rhs, uint_t level, DoFType flag = All);

  inline void prolongate_impl(uint_t sourceLevel, DoFType flag = All);

  inline void prolongateQuadratic_impl(uint_t sourceLevel, DoFType flag = All);

  inline void restrict_impl(uint_t sourceLevel, DoFType flag = All);

  inline void enumerate_impl(uint_t level, uint_t& num);
};

template< typename ValueType >
inline void P2Function< ValueType >::interpolate_impl(std::function< ValueType(const hhg::Point3D&) > & expr, uint_t level, DoFType flag)
{
  vertexDoFFunction.interpolate(expr, level, flag);
  edgeDoFFunction.interpolate(expr, level, flag);
}

template< typename ValueType >
inline void P2Function< ValueType >::assign_impl(const std::vector<ValueType> scalars, const std::vector<P2Function< ValueType >*> functions, size_t level, DoFType flag)
{
  std::vector<P1Function< ValueType >*> vertexFunctions;
  std::vector<EdgeDoFFunction< ValueType >*> edgeFunctions;

  for (auto function : functions) {
    vertexFunctions.push_back(&function->getVertexDoFFunction());
    edgeFunctions.push_back(&function->getEdgeDoFFunction());
  }

  vertexDoFFunction.assign(scalars, vertexFunctions, level, flag);
  edgeDoFFunction.assign(scalars, edgeFunctions, level, flag);
}

template< typename ValueType >
inline void P2Function< ValueType >::add_impl(const std::vector<ValueType> scalars, const std::vector<P2Function< ValueType >*> functions, size_t level, DoFType flag)
{
  std::vector<P1Function< ValueType >*> vertexFunctions;
  std::vector<EdgeDoFFunction< ValueType >*> edgeFunctions;

  for (auto function : functions) {
    vertexFunctions.push_back(&function->getVertexDoFFunction());
    edgeFunctions.push_back(&function->getEdgeDoFFunction());
  }

  vertexDoFFunction.add(scalars, vertexFunctions, level, flag);
  edgeDoFFunction.add(scalars, edgeFunctions, level, flag);
}

template< typename ValueType >
inline real_t P2Function< ValueType >::dot_impl(P2Function< ValueType >& rhs, size_t level, DoFType flag)
{
  real_t sum = vertexDoFFunction.dot(rhs.getVertexDoFFunction(), level, flag);
  sum += edgeDoFFunction.dot(rhs.getEdgeDoFFunction(), level, flag);
  return sum;
}

template< typename ValueType >
inline void P2Function< ValueType >::prolongate_impl(size_t sourceLevel, DoFType flag)
{
  WALBERLA_ABORT("not implemented");
}

template< typename ValueType >
inline void P2Function< ValueType >::prolongateQuadratic_impl(size_t sourceLevel, DoFType flag)
{
  WALBERLA_ABORT("not implemented");
}

template< typename ValueType >
inline void P2Function< ValueType >::restrict_impl(size_t sourceLevel, DoFType flag)
{
  WALBERLA_ABORT("not implemented");
}

template< typename ValueType >
inline void P2Function< ValueType >::enumerate_impl(uint_t level, uint_t& num)
{
  vertexDoFFunction.enumerate(level, num);
  edgeDoFFunction.enumerate(level, num);
}

}
