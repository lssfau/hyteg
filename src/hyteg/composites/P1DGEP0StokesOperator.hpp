/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"

namespace hyteg {
    namespace dg {
        namespace eg {
            template<typename VelocityBlockOperator, typename EnergyNormOperator>
            class EGP0StokesOperatorType : public Operator<EGP0StokesFunction<real_t>, EGP0StokesFunction<real_t> > {
            public:
                typedef VelocityBlockOperator VelocityBlockOperator_T;
                typedef EnergyNormOperator EnergyNormOperator_T;

                EGP0StokesOperatorType(const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
                                       size_t maxLevel)
                        : Operator(storage, minLevel, maxLevel), velocityBlockOp(storage, minLevel, maxLevel),
                          div(storage, minLevel, maxLevel), divT(storage, minLevel, maxLevel),
                          energyNormOp(storage, minLevel, maxLevel) {}

                void apply(const EGP0StokesFunction<real_t> &src,
                           const EGP0StokesFunction<real_t> &dst,
                           const uint_t level,
                           const DoFType flag) const {
                    velocityBlockOp.apply(src.uvw(), dst.uvw(), level, flag, Replace);
                    divT.apply(src.p(), dst.uvw(), level, flag, Add);
                    div.apply(src.uvw(), dst.p(), level, flag, Replace);
                }

                void toMatrix(const std::shared_ptr<SparseMatrixProxy> &mat,
                              const EGP0StokesFunction<idx_t> &src,
                              const EGP0StokesFunction<idx_t> &dst,
                              size_t level,
                              DoFType flag) const {
                    velocityBlockOp.toMatrix(mat, src.uvw(), dst.uvw(), level, flag);
                    divT.toMatrix(mat, src.p(), dst.uvw(), level, flag);
                    div.toMatrix(mat, src.uvw(), dst.p(), level, flag);
                }

                VelocityBlockOperator_T velocityBlockOp;
                EGToP0DivOperator div;
                P0ToEGDivTOperator divT;
                EnergyNormOperator_T energyNormOp;
            };

            typedef EGP0StokesOperatorType<EGSIPGLaplaceOperator, EGLaplaceEnergyNormOperator> EGP0StokesOperator;
            typedef EGP0StokesOperatorType<EGIIPGLaplaceOperator, EGLaplaceEnergyNormOperator> EGP0IIPGStokesOperator;


            template<typename VelocityBlockOperator>
            class EGP0EpsilonStokesOperatorType
                    : public Operator<EGP0StokesFunction<real_t>, EGP0StokesFunction<real_t> > {
            public:

                typedef VelocityBlockOperator VelocityBlockOperator_T;
                typedef EGEpsilonEnergyNormOperator EnergyNormOperator_T;

                EGP0EpsilonStokesOperatorType(const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel,
                                              size_t maxLevel, std::function<real_t(const Point3D &)> viscosity)
                        : Operator(storage, minLevel, maxLevel), eps(storage, minLevel, maxLevel, viscosity),
                          div(storage, minLevel, maxLevel), divT(storage, minLevel, maxLevel),
                          energyNormOp(storage, minLevel, maxLevel, viscosity) {}

                void apply(const EGP0StokesFunction<real_t> &src,
                           const EGP0StokesFunction<real_t> &dst,
                           const uint_t level,
                           const DoFType flag) const {
                    eps.apply(src.uvw(), dst.uvw(), level, flag, Replace);
                    divT.apply(src.p(), dst.uvw(), level, flag, Add);
                    div.apply(src.uvw(), dst.p(), level, flag, Replace);
                }

                void toMatrix(const std::shared_ptr<SparseMatrixProxy> &mat,
                              const EGP0StokesFunction<idx_t> &src,
                              const EGP0StokesFunction<idx_t> &dst,
                              size_t level,
                              DoFType flag) const {
                    eps.toMatrix(mat, src.uvw(), dst.uvw(), level, flag);
                    divT.toMatrix(mat, src.p(), dst.uvw(), level, flag);
                    div.toMatrix(mat, src.uvw(), dst.p(), level, flag);
                }

                VelocityBlockOperator_T eps;
                EGToP0DivOperator div;
                P0ToEGDivTOperator divT;
                EnergyNormOperator_T energyNormOp;
            };


            typedef EGP0EpsilonStokesOperatorType<EGEpsilonOperator> EGP0EpsilonStokesOperator;
        } // namespace eg
    } // namespace dg
} // namespace hyteg
