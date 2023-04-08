/*
* Copyright (c) 2017-2022 Nils Kohl.
*
* This file is part of HyTeG
* (see https://i10git.cs.fau.de/hyteg/hyteg).
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "hyteg/communication/Syncing.hpp"
#include <hyteg/communication/Syncing.hpp>
#include <hyteg/dgfunctionspace/DGVectorLaplaceForm.hpp>
#include <hyteg/egfunctionspace/EGEpsilonEnergyNormForm.hpp>

#include "hyteg/dgfunctionspace/DGDivForm.hpp"
#include "hyteg/dgfunctionspace/DGFormAbort.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGVectorMassForm.hpp"
#include "hyteg/dgfunctionspace/P0_to_P1_divt_form.hpp"
#include "hyteg/egfunctionspace/EGConstEpsilonForm.hpp"
#include "hyteg/egfunctionspace/EGDivForm.hpp"
#include "hyteg/egfunctionspace/EGEpsilonForm.hpp"
#include "hyteg/egfunctionspace/EGMassForm.hpp"

#include "hyteg/egfunctionspace/EGIIPGVectorLaplaceForm.hpp"
#include "hyteg/egfunctionspace/EGVectorLaplaceForm.hpp"
#include "hyteg/egfunctionspace/EGNIPGVectorLaplaceForm.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

namespace hyteg {

    using namespace dg;
    using facedof::FaceType;
    using indexing::Index;
    using volumedofspace::indexing::VolumeDoFMemoryLayout;
    using walberla::int_c;
    using walberla::real_t;

    template<typename Form>
    class P0ToP1Operator : public Operator<P0Function<real_t>, P1Function<real_t> > {
    public:
        P0ToP1Operator(const std::shared_ptr<PrimitiveStorage> &storage,
                       uint_t minLevel,
                       uint_t maxLevel,
                       std::shared_ptr<Form> form = std::make_shared<Form>())
                : Operator<P0Function<real_t>, P1Function<real_t> >(storage, minLevel, maxLevel), form_(form) {}

        typedef Form FormType;

        void apply(const P0Function<real_t> &src,
                   const P1Function<real_t> &dst,
                   size_t level,
                   DoFType flag,
                   UpdateType updateType) const override {
            assembleAndOrApply(src, dst, level, flag, nullptr, updateType);
        }

        void toMatrix(const std::shared_ptr<SparseMatrixProxy> &mat,
                      const P0Function<idx_t> &src,
                      const P1Function<idx_t> &dst,
                      size_t level,
                      DoFType flag) const override {
            assembleAndOrApply(src, dst, level, flag, mat, Replace);
        }

        template<typename VType>
        VType *p1Data(const P1Function<VType> &function,
                      const std::shared_ptr<PrimitiveStorage> &storage,
                      const PrimitiveID &pid,
                      uint_t level) const {
            if (storage->hasGlobalCells()) {
                WALBERLA_ASSERT(storage->cellExistsLocally(pid));
                auto cell = storage->getCell(pid);
                return cell->getData(function.getCellDataID())->getPointer(level);
            } else {
                WALBERLA_ASSERT(storage->faceExistsLocally(pid));
                auto face = storage->getFace(pid);
                return face->getData(function.getFaceDataID())->getPointer(level);
            }
        }

    private:
        /// \brief This is similar to the implementation in the dg::DGOperator class.
        template<typename VType>
        inline void assembleAndOrApply(const P0Function<VType> &src,
                                       const P1Function<VType> &dst,
                                       size_t level,
                                       DoFType flag,
                                       const std::shared_ptr<SparseMatrixProxy> &mat,
                                       UpdateType updateType = Replace) const {
            // To avoid code duplication in this already long method, the implementation "fuses" the 2D and 3D implementation.
            // This more or less serves as a reference - for better performance the matrix-vector multiplication should be specialized.

            DGBasisLinearLagrange_Example dstBasis;

            using indexing::Index;
            using volumedofspace::indexing::ElementNeighborInfo;

            // WALBERLA_CHECK( updateType == Replace );

            const auto storage = this->getStorage();

            int dim = 2;
            if (storage->hasGlobalCells()) {
                dim = 3;
            }

            std::vector<PrimitiveID> pids;
            if (dim == 2) {
                pids = storage->getFaceIDs();
            } else {
                pids = storage->getCellIDs();
            }

            src.communicate(level);

            if (updateType == Replace && mat == nullptr) {
                // We need to zero the destination array (including halos).
                // However, we must not zero out anything that is not flagged with the specified BCs.
                // Therefore we first zero out everything that flagged, and then, later,
                // the halos of the highest dim primitives.

                dst.interpolate(real_c(0), level, flag);
            }

            for (const auto &pid: pids) {
                const auto srcPolyDegree = 0;
                const auto dstPolyDegree = 1;

                const auto numSrcDofs = 1;
                const auto numDstDofs = dim + 1;

                const auto srcDofMemory = src.getDGFunction()->volumeDoFFunction()->dofMemory(pid, level);
                auto dstDofMemory = p1Data<VType>(dst, storage, pid, level);

                const auto srcMemLayout = src.getDGFunction()->volumeDoFFunction()->memoryLayout();

                std::map<uint_t, VType *> glMemory;

                if (dim == 2) {
                    WALBERLA_ASSERT(storage->faceExistsLocally(pid));
                    const auto face = storage->getFace(pid);
                    for (const auto &[n, _]: face->getIndirectNeighborFaceIDsOverEdges()) {
                        glMemory[n] = src.getDGFunction()->volumeDoFFunction()->glMemory(pid, level, n);
                    }
                } else {
                    WALBERLA_ASSERT(storage->cellExistsLocally(pid));
                    const auto cell = storage->getCell(pid);
                    for (const auto &[n, _]: cell->getIndirectNeighborCellIDsOverFaces()) {
                        glMemory[n] = src.getDGFunction()->volumeDoFFunction()->glMemory(pid, level, n);
                    }
                }

                // zero out halos for matrix-free application
                if (mat == nullptr) {
                    if (dim == 2) {
                        for (const auto &idx: vertexdof::macroface::Iterator(level)) {
                            if (vertexdof::macroface::isVertexOnBoundary(level, idx)) {
                                auto arrayIdx = vertexdof::macroface::index(level, idx.x(), idx.y());
                                dstDofMemory[arrayIdx] = real_c(0);
                            }
                        }
                    } else {
                        for (const auto &idx: vertexdof::macrocell::Iterator(level)) {
                            if (!vertexdof::macrocell::isOnCellFace(idx, level).empty()) {
                                auto arrayIdx = vertexdof::macrocell::index(level, idx.x(), idx.y(), idx.z());
                                dstDofMemory[arrayIdx] = real_c(0);
                            }
                        }
                    }
                }

                const uint_t numMicroVolTypes = (storage->hasGlobalCells() ? 6 : 2);

                for (uint_t microVolType = 0; microVolType < numMicroVolTypes; microVolType++) {
                    if (dim == 2 && microVolType >= 2) {
                        break;
                    }

                    auto faceType = facedof::allFaceTypes[microVolType];
                    auto cellType = celldof::allCellTypes[microVolType];

                    auto itFace = facedof::macroface::Iterator(level, faceType).begin();
                    auto itCell = celldof::macrocell::Iterator(level, cellType).begin();

                    while ((dim == 2 && itFace != itFace.end()) || (dim == 3 && itCell != itCell.end())) {
                        Index elementIdx;

                        if (dim == 2) {
                            elementIdx = *itFace;
                            itFace++;
                        } else {
                            elementIdx = *itCell;
                            itCell++;
                        }

                        // TODO: all these coord computations can be executed _once_ and then the coordinates can be incremented by h
                        // TODO: blending

                        // This object does the heavy lifting of computing all required coordinates and normals.
                        ElementNeighborInfo neighborInfo;

                        if (dim == 2) {
                            neighborInfo = ElementNeighborInfo(elementIdx, faceType, level, src.getBoundaryCondition(),
                                                               pid, storage_);
                        } else {
                            neighborInfo = ElementNeighborInfo(elementIdx, cellType, level, src.getBoundaryCondition(),
                                                               pid, storage_);
                        }

                        // We only write to the DoFs in the current volume, let's prepare a temporary vector for that.
                        Eigen::Matrix<real_t, Eigen::Dynamic, 1> dstDofs;
                        dstDofs.resize(numDstDofs, Eigen::NoChange_t::NoChange);
                        dstDofs.setZero();

                        /////////////////////////
                        // Volume contribution //
                        /////////////////////////

                        Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> localMat;
                        localMat.resize(numDstDofs, numSrcDofs);

                        // Little difference here is that the dst is now a CG P1 function.
                        // So we need to write the DoFs a little differently and set the basis manually.

                        form_->integrateVolume(dim,
                                               neighborInfo.elementVertexCoords(),
                                               *src.getDGFunction()->basis(),
                                               dstBasis,
                                               srcPolyDegree,
                                               dstPolyDegree,
                                               localMat);

                        // P0 has only one DoF
                        auto srcDoF = dim == 2 ? srcDofMemory[volumedofspace::indexing::index(
                                elementIdx.x(), elementIdx.y(), faceType, 0, 1, level, srcMemLayout)]
                                               : srcDofMemory[volumedofspace::indexing::index(
                                        elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, 0, 1, level,
                                        srcMemLayout)];


                        // Getting the vertex DoF indices for the current micro volume.
                        std::vector<Index> vertexDoFIndices;
                        if (dim == 2) {
                            auto vertexDoFIndicesArray = facedof::macroface::getMicroVerticesFromMicroFace(elementIdx,
                                                                                                           faceType);
                            vertexDoFIndices.insert(vertexDoFIndices.begin(), vertexDoFIndicesArray.begin(),
                                                    vertexDoFIndicesArray.end());
                        } else {
                            auto vertexDoFIndicesArray = celldof::macrocell::getMicroVerticesFromMicroCell(elementIdx,
                                                                                                           cellType);
                            vertexDoFIndices.insert(vertexDoFIndices.begin(), vertexDoFIndicesArray.begin(),
                                                    vertexDoFIndicesArray.end());
                        }

                        if (mat == nullptr) {
                            // Matrix-vector multiplication.
                            dstDofs += localMat * srcDoF;
                        } else {
                            // Sparse assembly.
                            for (uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++) {
                                if (dim == 2) {
                                    const auto globalRowIdx = dstDofMemory[vertexdof::macroface::index(
                                            level, vertexDoFIndices[dstDofIdx].x(), vertexDoFIndices[dstDofIdx].y())];
                                    const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                                            elementIdx.x(), elementIdx.y(), faceType, 0, 1, level, srcMemLayout)];
                                    mat->addValue(globalRowIdx, globalColIdx, localMat(Eigen::Index(dstDofIdx), 0));
                                } else {
                                    const auto globalRowIdx = dstDofMemory[vertexdof::macrocell::index(level,
                                                                                                       vertexDoFIndices[dstDofIdx].x(),
                                                                                                       vertexDoFIndices[dstDofIdx].y(),
                                                                                                       vertexDoFIndices[dstDofIdx].z())];
                                    const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                                            elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, 0, 1, level,
                                            srcMemLayout)];
                                    mat->addValue(globalRowIdx, globalColIdx, localMat(Eigen::Index(dstDofIdx), 0));
                                }
                            }
                        }

                        if (!form_->onlyVolumeIntegrals()) {
                            /////////////////////////////
                            // Interface contributions //
                            /////////////////////////////

                            // Loop over neighboring volumes.
                            for (uint_t n = 0; n < uint_c(dim + 1); n++) {
                                /////////////////////
                                // Domain boundary //
                                /////////////////////

                                if (neighborInfo.atMacroBoundary(n) &&
                                    neighborInfo.neighborBoundaryType(n) == DirichletBoundary) {
                                    ////////////////////////
                                    // Dirichlet boundary //
                                    ////////////////////////

                                    localMat.setZero();
                                    form_->integrateFacetDirichletBoundary(dim,
                                                                           neighborInfo.elementVertexCoords(),
                                                                           neighborInfo.interfaceVertexCoords(n),
                                                                           neighborInfo.oppositeVertexCoords(n),
                                                                           neighborInfo.outwardNormal(n),
                                                                           *src.getDGFunction()->basis(),
                                                                           dstBasis,
                                                                           srcPolyDegree,
                                                                           dstPolyDegree,
                                                                           localMat);

                                    if (mat == nullptr) {
                                        // Matrix-vector multiplication.
                                        dstDofs += localMat * srcDoF;
                                    } else {
                                        // Sparse assembly.
                                        for (uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++) {
                                            if (dim == 2) {
                                                const auto globalRowIdx = dstDofMemory[vertexdof::macroface::index(
                                                        level, vertexDoFIndices[dstDofIdx].x(),
                                                        vertexDoFIndices[dstDofIdx].y())];
                                                const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                                                        elementIdx.x(), elementIdx.y(), faceType, 0, 1, level,
                                                        srcMemLayout)];
                                                mat->addValue(globalRowIdx, globalColIdx,
                                                              localMat(Eigen::Index(dstDofIdx), 0));
                                            } else {
                                                const auto globalRowIdx =
                                                        dstDofMemory[vertexdof::macrocell::index(level,
                                                                                                 vertexDoFIndices[dstDofIdx].x(),
                                                                                                 vertexDoFIndices[dstDofIdx].y(),
                                                                                                 vertexDoFIndices[dstDofIdx].z())];
                                                const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                                                        elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, 0, 1,
                                                        level, srcMemLayout)];
                                                mat->addValue(globalRowIdx, globalColIdx,
                                                              localMat(Eigen::Index(dstDofIdx), 0));
                                            }
                                        }
                                    }
                                } else if (neighborInfo.atMacroBoundary(n) &&
                                           neighborInfo.neighborBoundaryType(n) == NeumannBoundary) {
                                    WALBERLA_ABORT("Neumann boundary handling not implemented.");
                                } else if (neighborInfo.atMacroBoundary(n) &&
                                           neighborInfo.neighborBoundaryType(n) == FreeslipBoundary) {
                                    WALBERLA_ABORT("Free-slip boundary handling not implemented.");
                                }

                                    //////////////////
                                    // Inner domain //
                                    //////////////////

                                else {
                                    ///////////////////////////////////
                                    // a) inner element contribution //
                                    ///////////////////////////////////

                                    localMat.setZero();
                                    form_->integrateFacetInner(dim,
                                                               neighborInfo.elementVertexCoords(),
                                                               neighborInfo.interfaceVertexCoords(n),
                                                               neighborInfo.oppositeVertexCoords(n),
                                                               neighborInfo.outwardNormal(n),
                                                               *src.getDGFunction()->basis(),
                                                               dstBasis,
                                                               srcPolyDegree,
                                                               dstPolyDegree,
                                                               localMat);

                                    if (mat == nullptr) {
                                        // Matrix-vector multiplication.
                                        dstDofs += localMat * srcDoF;
                                    } else {
                                        // Sparse assembly.
                                        for (uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++) {
                                            if (dim == 2) {
                                                const auto globalRowIdx = dstDofMemory[vertexdof::macroface::index(
                                                        level, vertexDoFIndices[dstDofIdx].x(),
                                                        vertexDoFIndices[dstDofIdx].y())];
                                                const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                                                        elementIdx.x(), elementIdx.y(), faceType, 0, 1, level,
                                                        srcMemLayout)];
                                                mat->addValue(globalRowIdx, globalColIdx,
                                                              localMat(Eigen::Index(dstDofIdx), 0));
                                            } else {
                                                const auto globalRowIdx =
                                                        dstDofMemory[vertexdof::macrocell::index(level,
                                                                                                 vertexDoFIndices[dstDofIdx].x(),
                                                                                                 vertexDoFIndices[dstDofIdx].y(),
                                                                                                 vertexDoFIndices[dstDofIdx].z())];
                                                const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                                                        elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, 0, 1,
                                                        level, srcMemLayout)];
                                                mat->addValue(globalRowIdx, globalColIdx,
                                                              localMat(Eigen::Index(dstDofIdx), 0));
                                            }
                                        }
                                    }

                                    ////////////////////////////////////////
                                    // b) coupling to neighboring element //
                                    ////////////////////////////////////////

                                    if (neighborInfo.atMacroBoundary(n) &&
                                        neighborInfo.neighborBoundaryType(n) == Inner) {

                                        ////////////////////////////////////////////////
                                        // i) micro-interface on macro-macro-boundary //
                                        ////////////////////////////////////////////////


                                        neighborInfo = neighborInfo.updateForMacroBoundary(n);

                                        localMat.setZero();
                                        form_->integrateFacetCoupling(dim,
                                                                      neighborInfo.elementVertexCoords(),
                                                                      neighborInfo.neighborElementVertexCoords( n ),
                                                                      neighborInfo.interfaceVertexCoords(n),
                                                                      neighborInfo.oppositeVertexCoords(n),
                                                                      neighborInfo.neighborOppositeVertexCoords( n ),
                                                                      neighborInfo.outwardNormal(n),
                                                                      *src.getDGFunction()->basis(),
                                                                      dstBasis,
                                                                      srcPolyDegree,
                                                                      dstPolyDegree,
                                                                      localMat);

                                        // Now we need the DoFs from the neighboring element.
                                        // There is only one DoF (the single P0/DG DoF)

                                        real_t nSrcDoF;
                                        uint_t nSrcDoFArrIndex = std::numeric_limits<uint_t>::max();

                                        if (dim == 2) {
                                            nSrcDoFArrIndex =
                                                    volumedofspace::indexing::indexNeighborInGhostLayer(
                                                            neighborInfo.macroBoundaryID(n),
                                                            elementIdx.x(),
                                                            elementIdx.y(),
                                                            faceType,
                                                            0,
                                                            numSrcDofs,
                                                            level,
                                                            srcMemLayout);
                                        } else {
                                            nSrcDoFArrIndex =
                                                    volumedofspace::indexing::indexNeighborInGhostLayer(
                                                            neighborInfo.macroBoundaryID(n),
                                                            elementIdx.x(),
                                                            elementIdx.y(),
                                                            elementIdx.z(),
                                                            cellType,
                                                            0,
                                                            numSrcDofs,
                                                            level,
                                                            srcMemLayout);
                                        }

                                        nSrcDoF = glMemory[neighborInfo.macroBoundaryID(n)][nSrcDoFArrIndex];

                                        if (mat == nullptr) {
                                            // Matrix-vector multiplication.
                                            dstDofs += localMat * nSrcDoF;
                                        } else {
                                            // Sparse assembly.
                                            // TODO: maybe there is a nicer way to do the gl stuff ...
                                            for (uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++) {
                                                uint_t globalRowIdx;
                                                if (dim == 2) {
                                                    globalRowIdx = dstDofMemory[vertexdof::macroface::index(
                                                            level, vertexDoFIndices[dstDofIdx].x(),
                                                            vertexDoFIndices[dstDofIdx].y())];
                                                } else {
                                                    globalRowIdx = dstDofMemory[vertexdof::macrocell::index(level,
                                                                                                            vertexDoFIndices[dstDofIdx].x(),
                                                                                                            vertexDoFIndices[dstDofIdx].y(),
                                                                                                            vertexDoFIndices[dstDofIdx].z())];
                                                }
                                                const auto globalColIdx = glMemory[neighborInfo.macroBoundaryID(
                                                        n)][nSrcDoFArrIndex];

                                                mat->addValue(globalRowIdx, globalColIdx, localMat(dstDofIdx, 0));
                                            }
                                        }
                                    } else {
                                        /////////////////////////////////////////
                                        // ii) micro-interface inside of macro //
                                        /////////////////////////////////////////

                                        localMat.setZero();
                                        form_->integrateFacetCoupling(dim,
                                                                      neighborInfo.elementVertexCoords(),
                                                                      neighborInfo.neighborElementVertexCoords(n),
                                                                      neighborInfo.interfaceVertexCoords(n),
                                                                      neighborInfo.oppositeVertexCoords(n),
                                                                      neighborInfo.neighborOppositeVertexCoords(n),
                                                                      neighborInfo.outwardNormal(n),
                                                                      *src.getDGFunction()->basis(),
                                                                      dstBasis,
                                                                      srcPolyDegree,
                                                                      dstPolyDegree,
                                                                      localMat);

                                        // Now we need the DoFs from the neighboring element.
                                        // P0 has only one DoF

                                        const auto nSrcDof =
                                                dim == 2 ?
                                                srcDofMemory[volumedofspace::indexing::index(
                                                        neighborInfo.neighborElementIndices(n).x(),
                                                        neighborInfo.neighborElementIndices(n).y(),
                                                        neighborInfo.neighborFaceType(n),
                                                        0,
                                                        1,
                                                        level,
                                                        srcMemLayout)] :
                                                srcDofMemory[volumedofspace::indexing::index(
                                                        neighborInfo.neighborElementIndices(n).x(),
                                                        neighborInfo.neighborElementIndices(n).y(),
                                                        neighborInfo.neighborElementIndices(n).z(),
                                                        neighborInfo.neighborCellType(n),
                                                        0,
                                                        1,
                                                        level,
                                                        srcMemLayout)];

                                        if (mat == nullptr) {
                                            // Matrix-vector multiplication.
                                            dstDofs += localMat * nSrcDof;
                                        } else {
                                            // TODO: improve this monster
                                            std::map<facedof::FaceType, uint_t> invFaceTypeMap;
                                            std::map<celldof::CellType, uint_t> invCellTypeMap;

                                            for (uint_t i = 0; i < 2; i++) {
                                                invFaceTypeMap[facedof::allFaceTypes[i]] = i;
                                            }
                                            for (uint_t i = 0; i < 6; i++) {
                                                invCellTypeMap[celldof::allCellTypes[i]] = i;
                                            }

                                            uint_t neighborMicroVolType;
                                            if (dim == 2) {
                                                neighborMicroVolType = invFaceTypeMap[neighborInfo.neighborFaceType(n)];
                                            } else {
                                                neighborMicroVolType = invCellTypeMap[neighborInfo.neighborCellType(n)];
                                            }

                                            // Sparse assembly.
                                            for (uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++) {
                                                if (dim == 2) {
                                                    const auto globalRowIdx = dstDofMemory[vertexdof::macroface::index(
                                                            level, vertexDoFIndices[dstDofIdx].x(),
                                                            vertexDoFIndices[dstDofIdx].y())];
                                                    const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                                                            neighborInfo.neighborElementIndices(n).x(),
                                                            neighborInfo.neighborElementIndices(n).y(),
                                                            neighborInfo.neighborFaceType(n),
                                                            0,
                                                            1,
                                                            level,
                                                            srcMemLayout)];
                                                    mat->addValue(globalRowIdx, globalColIdx,
                                                                  localMat(Eigen::Index(dstDofIdx), 0));
                                                } else {
                                                    const auto globalRowIdx =
                                                            dstDofMemory[vertexdof::macrocell::index(level,
                                                                                                     vertexDoFIndices[dstDofIdx].x(),
                                                                                                     vertexDoFIndices[dstDofIdx].y(),
                                                                                                     vertexDoFIndices[dstDofIdx].z())];
                                                    const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                                                            neighborInfo.neighborElementIndices(n).x(),
                                                            neighborInfo.neighborElementIndices(n).y(),
                                                            neighborInfo.neighborElementIndices(n).z(),
                                                            neighborInfo.neighborCellType(n),
                                                            0,
                                                            1,
                                                            level,
                                                            srcMemLayout)];
                                                    mat->addValue(globalRowIdx, globalColIdx,
                                                                  localMat(Eigen::Index(dstDofIdx), 0));
                                                }
                                            }
                                        }
                                    }
                                }
                            } // End loop over neighboring volumes.
                        }    // End if( !onlyVolumeIntegrals() )

                        if (mat == nullptr) {
                            // Write DoFs.
                            for (uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++) {
                                if (dim == 2) {
                                    dstDofMemory[vertexdof::macroface::index(
                                            level, vertexDoFIndices[dstDofIdx].x(),
                                            vertexDoFIndices[dstDofIdx].y())] += dstDofs(dstDofIdx);
                                } else {
                                    dstDofMemory[vertexdof::macrocell::index(level,
                                                                             vertexDoFIndices[dstDofIdx].x(),
                                                                             vertexDoFIndices[dstDofIdx].y(),
                                                                             vertexDoFIndices[dstDofIdx].z())] += dstDofs(
                                            dstDofIdx);
                                }
                            }
                        }
                    }
                }
            }

            if (mat == nullptr) {
                if (dim == 2) {
                    dst.template communicateAdditively<Face, Edge>(level, DoFType::All ^ flag, *storage_,
                                                                   updateType == Replace);
                    dst.template communicateAdditively<Face, Vertex>(level, DoFType::All ^ flag, *storage_,
                                                                     updateType == Replace);
                } else {
                    dst.template communicateAdditively<Cell, Face>(level, DoFType::All ^ flag, *storage_,
                                                                   updateType == Replace);
                    dst.template communicateAdditively<Cell, Edge>(level, DoFType::All ^ flag, *storage_,
                                                                   updateType == Replace);
                    dst.template communicateAdditively<Cell, Vertex>(level, DoFType::All ^ flag, *storage_,
                                                                     updateType == Replace);
                }
            }
        }

        std::shared_ptr<Form> form_;
    };

    // P1toP0 Stokes divergence
    typedef P0ToP1Operator<dg::p0_to_p1_divt_0_affine_q0> P0ToP1ConstantDivTxOperator;
    typedef P0ToP1Operator<dg::p0_to_p1_divt_1_affine_q0> P0ToP1ConstantDivTyOperator;
    typedef P0ToP1Operator<dg::p0_to_p1_divt_2_affine_q0> P0ToP1ConstantDivTzOperator;

    // EG Laplace operator couplings with different DG schemes
    typedef P0ToP1Operator<dg::eg::EGVectorLaplaceForm_P1E_0> EGVectorLaplaceP0ToP1Coupling_X;
    typedef P0ToP1Operator<dg::eg::EGVectorLaplaceForm_P1E_1> EGVectorLaplaceP0ToP1Coupling_Y;
    typedef P0ToP1Operator<dg::eg::EGVectorLaplaceForm_P1E_2> EGVectorLaplaceP0ToP1Coupling_Z;


    typedef P0ToP1Operator<dg::eg::EGNIPGVectorLaplaceFormP1E_0> EGNIPGVectorLaplaceP0ToP1Coupling_X;
    typedef P0ToP1Operator<dg::eg::EGNIPGVectorLaplaceFormP1E_1> EGNIPGVectorLaplaceP0ToP1Coupling_Y;
    typedef P0ToP1Operator<dg::eg::EGNIPGVectorLaplaceFormP1E_2> EGNIPGVectorLaplaceP0ToP1Coupling_Z;

    typedef P0ToP1Operator<dg::eg::EGIIPGVectorLaplaceFormP1E_0> EGIIPGVectorLaplaceP0ToP1Coupling_X;
    typedef P0ToP1Operator<dg::eg::EGIIPGVectorLaplaceFormP1E_1> EGIIPGVectorLaplaceP0ToP1Coupling_Y;
    typedef P0ToP1Operator<dg::eg::EGIIPGVectorLaplaceFormP1E_2> EGIIPGVectorLaplaceP0ToP1Coupling_Z;

    // EG Epsilon operator couplings
    typedef P0ToP1Operator<dg::eg::EGConstEpsilonFormP1E_0> EGConstantEpsilonP0ToP1Coupling_X;
    typedef P0ToP1Operator<dg::eg::EGConstEpsilonFormP1E_1> EGConstantEpsilonP0ToP1Coupling_Y;
    typedef P0ToP1Operator<dg::eg::EGConstEpsilonFormP1E_2> EGConstantEpsilonP0ToP1Coupling_Z;

    typedef P0ToP1Operator<dg::eg::EGEpsilonEnergyNormFormP1E_0> EGEpsilonEnergyNormP0ToP1Coupling_X;
    typedef P0ToP1Operator<dg::eg::EGEpsilonEnergyNormFormP1E_1> EGEpsilonEnergyNormP0ToP1Coupling_Y;
    typedef P0ToP1Operator<dg::eg::EGEpsilonEnergyNormFormP1E_2> EGEpsilonEnergyNormP0ToP1Coupling_Z;

    typedef P0ToP1Operator<dg::eg::EGEpsilonFormP1E_0> EGEpsilonP0ToP1Coupling_X;
    typedef P0ToP1Operator<dg::eg::EGEpsilonFormP1E_1> EGEpsilonP0ToP1Coupling_Y;
    typedef P0ToP1Operator<dg::eg::EGEpsilonFormP1E_2> EGEpsilonP0ToP1Coupling_Z;

// EG mass couplings
    typedef P0ToP1Operator<dg::eg::EGVectorMassFormP1E_0> EGMassP0ToP1Coupling_X;
    typedef P0ToP1Operator<dg::eg::EGVectorMassFormP1E_1> EGMassP0ToP1Coupling_Y;
    typedef P0ToP1Operator<dg::eg::EGVectorMassFormP1E_2> EGMassP0ToP1Coupling_Z;
//typedef P0ToP1Operator< dg::DGFormAbort> EGMassP1ToP0Coupling_Z;

//typedef P0ToP1Operator< dg::eg::EGDivFormP1E > P0ToP1ConstantP1EDGDivergenceCouplingOperator;

} // namespace hyteg