#pragma once

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

namespace hhg {

class EdgeDoFToVertexDoFOperator : public Operator<hhg::EdgeDoFFunction< real_t >, hhg::P1Function < real_t > >
{
public:

EdgeDoFToVertexDoFOperator(const std::shared_ptr <PrimitiveStorage> &storage, size_t minLevel, size_t maxLevel);

void apply_impl(EdgeDoFFunction< real_t >& src,  P1Function< real_t >& dst, uint_t level, DoFType flag, UpdateType updateType) final;

};

}// namespace hhg