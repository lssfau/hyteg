
#pragma once

#include <fmt/format.h>
#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "P1DataHandling.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics.hpp"

#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_pspg.h"
#include "tinyhhg_core/p1functionspace/generated/p1_stokes_epsilon.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "tinyhhg_core/p1functionspace/P1Memory.hpp"

#include "P1Vertex.hpp"
#include "P1Edge.hpp"
#include "P1Face.hpp"

namespace hhg
{

template<class UFCOperator,  bool Diagonal = false>
class P1ElementwiseOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  P1ElementwiseOperator(const std::shared_ptr< PrimitiveStorage > & storage, const std::array<const hhg::P1Function<real_t>*, 2>& coords, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel), coords_(coords)
  {
    computeElementMatrix_ = [this](Matrix3r& matrix, const real_t coords[6]) { this->computeElementMatrix(matrix, coords); };
  }

  ~P1ElementwiseOperator()
  {
  }

private:


  void apply_impl(P1Function< real_t >& src, P1Function< real_t >& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    std::array<const PrimitiveDataID<FunctionMemory<real_t>, Vertex>, 2> vertexCoordIds{{coords_[0]->getVertexDataID(), coords_[1]->getVertexDataID()}};
    std::array<const PrimitiveDataID<FunctionMemory<real_t>, Edge>, 2> edgeCoordIds{{coords_[0]->getEdgeDataID(), coords_[1]->getEdgeDataID()}};
    std::array<const PrimitiveDataID<FunctionMemory<real_t>, Face>, 2> faceCoordIds{{coords_[0]->getFaceDataID(), coords_[1]->getFaceDataID()}};

    src.getCommunicator(level)->startCommunication<Edge, Vertex>();
    src.getCommunicator(level)->startCommunication<Face, Edge>();
    src.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        vertexdof::macrovertex::applyElementwise< real_t >(level, vertex, storage_, computeElementMatrix_, src.getVertexDataID(), dst.getVertexDataID(), vertexCoordIds, updateType);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    src.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        P1Edge::applyElementwise< real_t >(level, edge, storage_, computeElementMatrix_, src.getEdgeDataID(), dst.getEdgeDataID(), edgeCoordIds, updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        vertexdof::macroface::applyElementwise< real_t >(level, face, computeElementMatrix_, src.getFaceDataID(), dst.getFaceDataID(), faceCoordIds, updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_gs_impl(P1Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("To be implemented")
  }

  void smooth_jac_impl(P1Function< real_t >& dst, P1Function< real_t >& rhs, P1Function< real_t >& tmp, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("To be implemented")
  }

#ifdef HHG_BUILD_WITH_PETSC
  void createMatrix_impl(P1Function< real_t >& src, P1Function< real_t >& dst, Mat& mat, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("To be implemented")
  }
#endif

private:

  void computeElementMatrix(Matrix3r& matrix, const real_t coords[6]) {
    ufcOperator_.tabulate_tensor(matrix.data(), NULL, coords, 0);
  }

  std::array<const hhg::P1Function<real_t>*, 2> coords_;
  UFCOperator ufcOperator_;

  std::function<void(Matrix3r&, const real_t[6])> computeElementMatrix_;
};

typedef P1ElementwiseOperator<p1_diffusion_cell_integral_0_otherwise> P1ElementwiseLaplaceOperator;
typedef P1ElementwiseOperator<p1_mass_cell_integral_0_otherwise> P1ElementwiseMassOperator;

}


