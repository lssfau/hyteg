#pragma once

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2Elements.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/forms/form_fenics_base/P2FenicsForm.hpp"
#include "tinyhhg_core/forms/form_fenics_base/P2ToP1FenicsForm.hpp"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t;

template< class P2ToP1Form >
class P2ToP1ConstantOperator : public Operator<P2Function < real_t>, P1Function<real_t> > {
 public:

  P2ToP1ConstantOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel),
        vertexToVertex(storage, minLevel, maxLevel),
        edgeToVertex(storage, minLevel, maxLevel)
  {
  }

  P1ConstantOperator< P2ToP1Form > const & getVertexToVertexOpr() const {
    return vertexToVertex;
  }

  EdgeDoFToVertexDoFOperator< P2ToP1Form > const & getEdgeToVertexOpr() const {
    return edgeToVertex;
  }

   void apply( const P2Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      vertexToVertex.apply( src.getVertexDoFFunction(), dst, level, flag, updateType );
      edgeToVertex.apply( src.getEdgeDoFFunction(), dst, level, flag, Add );
   }

  void smooth_gs(P2Function< real_t > & dst, P1Function< real_t > & rhs, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("not implemented");
  }

private:

  P1ConstantOperator< P2ToP1Form > vertexToVertex;
  EdgeDoFToVertexDoFOperator< P2ToP1Form > edgeToVertex;

};

typedef P2ToP1ConstantOperator< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > > P2ToP1ConstantDivxOperator;
typedef P2ToP1ConstantOperator< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > > P2ToP1ConstantDivyOperator;
typedef P2ToP1ConstantOperator< P2ToP1FenicsForm< fenics::NoAssemble,                     p2_to_p1_tet_div_tet_cell_integral_2_otherwise > > P2ToP1ConstantDivzOperator;

}
