#pragma once

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_pspg.h"
#include "tinyhhg_core/p1functionspace/generated/p1_stokes_epsilon.h"
#include "tinyhhg_core/p1functionspace/generated/p1_polar_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_polar_laplacian.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include <array>

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"

#include "tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"

#include "tinyhhg_core/p1functionspace/elementwise/VertexDoFMacroVertexElementwise.hpp"
#include "tinyhhg_core/p1functionspace/elementwise/VertexDoFMacroEdgeElementwise.hpp"
#include "tinyhhg_core/p1functionspace/elementwise/VertexDoFMacroFaceElementwise.hpp"

#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "P1DataHandling.hpp"

namespace hhg {

  template < class UFCOperator, bool Diagonal = false >
  class P1ElementwiseOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
  {
  public:
    P1ElementwiseOperator( const std::shared_ptr< PrimitiveStorage >&               storage,
                           const std::array< const hhg::P1Function< real_t >*, 2 >& coords,
                           size_t                                                   minLevel,
                           size_t                                                   maxLevel )
      : Operator( storage, minLevel, maxLevel )
      , coords_( coords )
    {
      computeElementMatrix_ = [this]( Matrix3r& matrix, const real_t coords[6] ) {
        this->computeElementMatrix( matrix, coords );
      };
    }

    ~P1ElementwiseOperator() {}

  private:
    void apply_impl( P1Function< real_t >& src,
                     P1Function< real_t >& dst,
                     size_t                level,
                     DoFType               flag,
                     UpdateType            updateType = Replace )
    {
      std::array< const PrimitiveDataID< FunctionMemory< real_t >, Vertex >, 2 > vertexCoordIDs{
        {coords_[0]->getVertexDataID(), coords_[1]->getVertexDataID()}};
      std::array< const PrimitiveDataID< FunctionMemory< real_t >, Edge >, 2 > edgeCoordIDs{
        {coords_[0]->getEdgeDataID(), coords_[1]->getEdgeDataID()}};
      std::array< const PrimitiveDataID< FunctionMemory< real_t >, Face >, 2 > faceCoordIDs{
        {coords_[0]->getFaceDataID(), coords_[1]->getFaceDataID()}};

      src.communicate< Vertex, Edge>( level );
      src.communicate< Edge, Face>( level );
      src.communicate< Face, Edge>( level );
      src.communicate< Edge, Vertex>( level );


      for( auto& it : storage_->getVertices() )
        {
          Vertex& vertex = *it.second;

          const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
          if( testFlag( vertexBC, flag ) )
            {
              vertexdof::macrovertex::applyElementwise< real_t >( level,
                                                                  vertex,
                                                                  storage_,
                                                                  computeElementMatrix_,
                                                                  src.getVertexDataID(),
                                                                  dst.getVertexDataID(),
                                                                  vertexCoordIDs,
                                                                  updateType );
            }
        }



      for( auto& it : storage_->getEdges() )
        {
          Edge& edge = *it.second;

          const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
          if( testFlag( edgeBC, flag ) )
            {
              vertexdof::macroedge::applyElementwise< real_t >( level,
                                                                edge,
                                                                storage_,
                                                                computeElementMatrix_,
                                                                src.getEdgeDataID(),
                                                                dst.getEdgeDataID(),
                                                                edgeCoordIDs,
                                                                updateType );
            }
        }

      for( auto& it : storage_->getFaces() )
        {
          Face& face = *it.second;

          const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
          if( testFlag( faceBC, flag ) )
            {
              vertexdof::macroface::applyElementwise< real_t >( level, face,
                                                                computeElementMatrix_,
                                                                src.getFaceDataID(),
                                                                dst.getFaceDataID(),
                                                                faceCoordIDs, updateType );
            }
        }

    }

    void smooth_gs_impl( P1Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag )
    {
      smooth_sor_impl( dst, rhs, real_t(1.0), level, flag );
    }

    void smooth_sor_impl( P1Function< real_t >& dst, P1Function< real_t >& rhs, real_t relax, size_t level, DoFType flag )
    {

      // extract handles to primitive coordinates
      std::array< const PrimitiveDataID< FunctionMemory< real_t >, Vertex >, 2 > vertexCoordIDs{
        { coords_[0]->getVertexDataID(), coords_[1]->getVertexDataID() }};
      std::array< const PrimitiveDataID< FunctionMemory< real_t >, Edge >, 2 > edgeCoordIDs{
        { coords_[0]->getEdgeDataID(), coords_[1]->getEdgeDataID() }};
      std::array< const PrimitiveDataID< FunctionMemory< real_t >, Face >, 2 > faceCoordIDs{
        { coords_[0]->getFaceDataID(), coords_[1]->getFaceDataID() }};

      // start pulling vertex halos
      dst.startCommunication< Edge, Vertex >( level );

      // start pulling edge halos
      dst.startCommunication< Face, Edge >( level );

      // end pulling vertex halos
      dst.endCommunication< Edge, Vertex >( level );

      for( auto& it : storage_->getVertices() )
        {
          Vertex& vertex = *it.second;

          const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
          if( testFlag( vertexBC, flag ) )
            {
              vertexdof::macrovertex::point_smooth_elementwise< real_t >( level, vertex,
                                                                          storage_,
                                                                          computeElementMatrix_,
                                                                          dst.getVertexDataID(),
                                                                          rhs.getVertexDataID(),
                                                                          dst.getVertexDataID(),
                                                                          vertexCoordIDs,
                                                                          relax );
            }
        }

      dst.startCommunication< Vertex, Edge >( level );

      // end pulling edge halos
      dst.endCommunication< Face, Edge >( level );

      for( auto& it : storage_->getEdges() )
        {
          Edge& edge = *it.second;

          const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
          if( testFlag( edgeBC, flag ) )
            {
              vertexdof::macroedge::point_smooth_elementwise< real_t >( level, edge,
                                                                        storage_,
                                                                        computeElementMatrix_,
                                                                        dst.getEdgeDataID(),
                                                                        rhs.getEdgeDataID(),
                                                                        dst.getEdgeDataID(),
                                                                        edgeCoordIDs,
                                                                        relax );
            }
        }

      dst.endCommunication< Vertex, Edge >( level );

      dst.startCommunication< Edge, Face >( level );

      for( auto& it : storage_->getFaces() )
        {
          Face& face = *it.second;

          const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
          if( testFlag( faceBC, flag ) )
            {
              vertexdof::macroface::point_smooth_elementwise< real_t >( level, face,
                                                                        computeElementMatrix_,
                                                                        dst.getFaceDataID(),
                                                                        rhs.getFaceDataID(),
                                                                        dst.getFaceDataID(),
                                                                        faceCoordIDs, relax );
            }
        }

      dst.endCommunication< Edge, Face >( level );
    }


    void smooth_jac_impl( P1Function< real_t >& dst,
                          P1Function< real_t >& rhs,
                          P1Function< real_t >& src,
                          size_t                level,
                          DoFType               flag )
    {

      // no weighting currently
      real_t relax(1.0);

      // extract handles to primitive coordinates
      std::array< const PrimitiveDataID< FunctionMemory< real_t >, Vertex >, 2 > vertexCoordIDs{
        { coords_[0]->getVertexDataID(), coords_[1]->getVertexDataID() }};
      std::array< const PrimitiveDataID< FunctionMemory< real_t >, Edge >, 2 > edgeCoordIDs{
        { coords_[0]->getEdgeDataID(), coords_[1]->getEdgeDataID() }};
      std::array< const PrimitiveDataID< FunctionMemory< real_t >, Face >, 2 > faceCoordIDs{
        { coords_[0]->getFaceDataID(), coords_[1]->getFaceDataID() }};

      // start pulling vertex halos
      src.startCommunication< Edge, Vertex >( level );

      // start pulling edge halos
      src.startCommunication< Face, Edge >( level );

      // end pulling vertex halos
      src.endCommunication< Edge, Vertex >( level );

      for( auto& it : storage_->getVertices() )
      {
        Vertex& vertex = *it.second;

        const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
        if( testFlag( vertexBC, flag ) )
          {
            vertexdof::macrovertex::point_smooth_elementwise< real_t >( level, vertex,
                                                                        storage_,
                                                                        computeElementMatrix_,
                                                                        dst.getVertexDataID(),
                                                                        rhs.getVertexDataID(),
                                                                        src.getVertexDataID(),
                                                                        vertexCoordIDs,
                                                                        relax );
          }
      }

      dst.startCommunication< Vertex, Edge >( level );

      // end pulling edge halos
      src.endCommunication< Face, Edge >( level );

      for( auto& it : storage_->getEdges() )
        {
          Edge& edge = *it.second;

          const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
          if( testFlag( edgeBC, flag ) )
            {
              vertexdof::macroedge::point_smooth_elementwise< real_t >( level, edge,
                                                                        storage_,
                                                                        computeElementMatrix_,
                                                                        dst.getEdgeDataID(),
                                                                        rhs.getEdgeDataID(),
                                                                        src.getEdgeDataID(),
                                                                        edgeCoordIDs,
                                                                        relax );
            }
        }

      dst.endCommunication< Vertex, Edge >( level );

      dst.startCommunication< Edge, Face >( level );

      for( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

        const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
           vertexdof::macroface::point_smooth_elementwise< real_t >( level, face,
                                                                     computeElementMatrix_,
                                                                     dst.getFaceDataID(),
                                                                     rhs.getFaceDataID(),
                                                                     src.getFaceDataID(),
                                                                     faceCoordIDs, relax );
         }
      }

      dst.endCommunication< Edge, Face >( level );
    }

#ifdef HHG_BUILD_WITH_PETSC
    void createMatrix_impl( P1Function< real_t >& src, P1Function< real_t >& dst, Mat& mat, size_t level, DoFType flag )
    {
      WALBERLA_ABORT( "To be implemented" )
        }
#endif

  private:
    void computeElementMatrix( Matrix3r& matrix, const real_t coords[6] )
    {
      ufcOperator_.tabulate_tensor( matrix.data(), NULL, coords, 0 );
    }

    std::array< const hhg::P1Function< real_t >*, 2 > coords_;
    UFCOperator                                       ufcOperator_;

    std::function< void( Matrix3r&, const real_t[6] ) > computeElementMatrix_;
  };

  typedef P1ElementwiseOperator< p1_diffusion_cell_integral_0_otherwise  > P1ElementwiseLaplaceOperator;
  typedef P1ElementwiseOperator< p1_mass_cell_integral_0_otherwise       > P1ElementwiseMassOperator;
  typedef P1ElementwiseOperator< p1_polar_mass_cell_integral_0_otherwise > P1ElementwisePolarMassOperator;
  typedef P1ElementwiseOperator< p1_polar_laplacian_cell_integral_0_otherwise > P1ElementwisePolarLaplaceOperator;

} // namespace hhg
