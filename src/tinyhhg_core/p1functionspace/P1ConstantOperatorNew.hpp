#pragma once

#include <array>

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"

#include "tinyhhg_core/p1functionspace/generated_new/P1FenicsForm.hpp"

namespace hhg {

using walberla::real_t;

template < class P1Form, bool Diagonal = false, bool Lumped = false, bool InvertDiagonal = false >
class P1ConstantOperatorNew : public Operator< P1Function< real_t >, P1Function< real_t > >
{
 public:
   P1ConstantOperatorNew( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );

   ~P1ConstantOperatorNew() override = default;

   void scale( real_t scalar );

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const;

   void smooth_gs( const P1Function< real_t >& dst, const P1Function< real_t >& rhs, size_t level, DoFType flag ) const;

   void smooth_sor( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag ) const;

   void smooth_jac( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    const P1Function< real_t >& tmp,
                    size_t                      level,
                    DoFType                     flag ) const;

   const PrimitiveDataID< StencilMemory< real_t >, Vertex >& getVertexStencilID() const { return vertexStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Edge >& getEdgeStencilID() const { return edgeStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Face >& getFaceStencilID() const { return faceStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Cell >& getCellStencilID() const { return cellStencilID_; }

 private:
   void assembleStencils();

   void assembleStencils3D();

 private:
   PrimitiveDataID< StencilMemory< real_t >, Vertex > vertexStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Edge >   edgeStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >   faceStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Cell >   cellStencilID_;

   P1Form form;
};

typedef P1ConstantOperatorNew< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise > > P1ConstantLaplaceOperatorNew;

} // namespace hhg
