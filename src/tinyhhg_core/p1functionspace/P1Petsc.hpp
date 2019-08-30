#pragma once

#include <tinyhhg_core/p1functionspace/P1Function.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFFunction.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp>

namespace hyteg {
namespace petsc {

inline void createVectorFromFunction(const P1Function<PetscScalar> &function,
                                     const P1Function<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  for (auto &it : function.getStorage()->getVertices()) {
    Vertex &vertex = *it.second;

    const DoFType vertexBC = function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if (testFlag(vertexBC, flag)) {
      vertexdof::macrovertex::createVectorFromFunction<PetscScalar>(vertex, function.getVertexDataID(), numerator.getVertexDataID(), vec, level);
    }
  }

  for (auto &it : function.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag)) {
      vertexdof::macroedge::createVectorFromFunction<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag)) {
      vertexdof::macroface::createVectorFromFunction<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getCells()) {
    Cell & cell = *it.second;

    const DoFType cellBC = function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if (testFlag(cellBC, flag)) {
      vertexdof::macrocell::createVectorFromFunction<PetscScalar>(level, cell, function.getCellDataID(), numerator.getCellDataID(), vec);
    }
  }
}

inline void createFunctionFromVector(const P1Function<PetscScalar> &function,
                                     const P1Function<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  for (auto &it : function.getStorage()->getVertices()) {
    Vertex &vertex = *it.second;

    const DoFType vertexBC = function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if (testFlag(vertexBC, flag)) {
      vertexdof::macrovertex::createFunctionFromVector<PetscScalar>(vertex, function.getVertexDataID(), numerator.getVertexDataID(), vec, level);
    }
  }

  function.startCommunication<Vertex, Edge>( level );
  function.endCommunication<Vertex, Edge>( level );

  for (auto &it : function.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag)) {
      vertexdof::macroedge::createFunctionFromVector<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  function.startCommunication<Edge, Face>( level );
  function.endCommunication<Edge, Face>( level );

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag)) {
      vertexdof::macroface::createFunctionFromVector<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }

  function.startCommunication<Face, Cell>( level );
  function.endCommunication<Face, Cell>( level );

  for (auto &it : function.getStorage()->getCells()) {
    Cell & cell = *it.second;

    const DoFType cellBC = function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if ( testFlag(cellBC, flag)) {
      vertexdof::macrocell::createFunctionFromVector<PetscScalar>(level, cell, function.getCellDataID(), numerator.getCellDataID(), vec);
    }
  }
}

inline void applyDirichletBC(const P1Function<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {
  for (auto &it : numerator.getStorage()->getVertices()) {
    Vertex &vertex = *it.second;

    const DoFType vertexBC = numerator.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if (testFlag(vertexBC, DirichletBoundary)) {
      vertexdof::macrovertex::applyDirichletBC(vertex, mat, level, numerator.getVertexDataID());
    }
  }

  for (auto &it : numerator.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = numerator.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, DirichletBoundary)) {
      vertexdof::macroedge::applyDirichletBC(level, edge, mat, numerator.getEdgeDataID());
    }
  }

  for (auto &it : numerator.getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = numerator.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, DirichletBoundary)) {
      vertexdof::macroface::applyDirichletBC(level, face, mat, numerator.getFaceDataID());
    }
  }

}

template < class OperatorType >
inline void createMatrix( const OperatorType&           opr,
                          const P1Function< PetscInt >& src,
                          const P1Function< PetscInt >& dst,
                          Mat&                          mat,
                          size_t                        level,
                          DoFType                       flag )
{
  const auto storage = src.getStorage();

  for (auto& it : opr.getStorage()->getVertices()) {
    Vertex& vertex = *it.second;

    const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if (testFlag(vertexBC, flag))
    {
      vertexdof::macrovertex::saveOperator(vertex, opr.getVertexStencilID(), src.getVertexDataID(), dst.getVertexDataID(), mat, level);
    }
  }

  for (auto& it : opr.getStorage()->getEdges()) {
    Edge& edge = *it.second;

    const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag))
    {
      vertexdof::macroedge::saveOperator(level, edge, *storage, opr.getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat);
    }
  }

  for (auto& it : opr.getStorage()->getFaces()) {
    Face& face = *it.second;

    const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag))
    {
      if ( storage->hasGlobalCells() )
      {
        vertexdof::macroface::saveOperator3D(level, face, *storage, opr.getFaceStencil3DID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
      }
      else
      {
        vertexdof::macroface::saveOperator(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
      }

    }
  }

  for (auto& it : opr.getStorage()->getCells()) {
    Cell & cell = *it.second;

    const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if (testFlag(cellBC, flag))
    {
      vertexdof::macrocell::saveOperator(level, cell, opr.getCellStencilID(), src.getCellDataID(), dst.getCellDataID(), mat);
    }
  }
}

}
}
