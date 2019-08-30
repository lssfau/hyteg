#pragma once

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2Elements.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/forms/form_fenics_base/P2FenicsForm.hpp"
#include "tinyhhg_core/forms/form_fenics_base/P1ToP2FenicsForm.hpp"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t;

template< class P1ToP2Form >
class P1ToP2ConstantOperator : public Operator<P1Function < real_t>, P2Function<real_t> > {
 public:

  P1ToP2ConstantOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel), vertexToVertex(storage, minLevel, maxLevel),
        vertexToEdge(storage, minLevel, maxLevel)
  {
  }


  P1ConstantOperator< P1ToP2Form > const & getVertexToVertexOpr() const {
    return vertexToVertex;
  }

  VertexDoFToEdgeDoFOperator< P1ToP2Form > const & getVertexToEdgeOpr() const {
    return vertexToEdge;
  }

   void apply( const P1Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      vertexToVertex.apply( src, dst.getVertexDoFFunction(), level, flag, updateType );
      vertexToEdge.apply( src, dst.getEdgeDoFFunction(), level, flag, updateType );
  }

  void smooth_gs(P1Function< real_t > & dst, P2Function< real_t > & rhs, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("not implemented");
  }

private:

  P1ConstantOperator< P1ToP2Form > vertexToVertex;
  VertexDoFToEdgeDoFOperator< P1ToP2Form > vertexToEdge;

};

typedef P1ToP2ConstantOperator< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > > P1ToP2ConstantDivTxOperator;
typedef P1ToP2ConstantOperator< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > > P1ToP2ConstantDivTyOperator;
typedef P1ToP2ConstantOperator< P1ToP2FenicsForm< fenics::NoAssemble,                      p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > > P1ToP2ConstantDivTzOperator;

}
