/*
 * Copyright (c) 2017-2021 Marcus Mohr.
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
#include <type_traits>

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScSparseMatrixProxy.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"
#include "mixed_operator/VectorLaplaceOperator.hpp"
namespace hyteg {

using walberla::real_t;

template < class VecFunctionType, class LaplOpType, class DivOpType, class DivOpTType > // different discretizations for Laplace operator and Divergence
class AugmentedLagrangianOperator : public Operator<  VecFunctionType,  VecFunctionType >
{

    public:
    
        AugmentedLagrangianOperator(
            const std::shared_ptr< PrimitiveStorage >& storage,
            size_t                                     level, // only one level, GKB not to be leveled any time soon
            real_t                                     Nu = 0 // nu = Augmented Laplacian Parameter
        ) 
        : Operator< VecFunctionType, VecFunctionType >( storage, level, level )
        , storage_(storage)
        , level_(level)
        , nu(Nu)
        , Lapl( storage, level, level )
        , div( storage, level, level )
        , divT( storage, level, level )
        { 
        }

        void apply( 
            const VecFunctionType&         src,
            const VecFunctionType&         dst,
            size_t                         level,
            DoFType                        flag,
            UpdateType                     updateType = Replace 
        ) const
        {
            // apply Augmented Lagrangian opterator:
            // ALOPx = (Laplace + nu Div^T Div)x
            // in [AR13] paper notation: M = W + nuAA^T
            if(nu > 0) {
                
            VecFunctionType tmp0("tmp0",storage_, level, level);
            typename DivOpType::dstType tmp1("tmp1",storage_, level, level);
            VecFunctionType tmp2("tmp2",storage_, level, level);
            Lapl.apply( src, tmp0, level, flag );
            div.apply( src, tmp1, level, flag );
            divT.apply( tmp1, tmp2, level, flag );
            dst.assign( {1, nu}, {tmp0, tmp2}, level,flag );
            } else {
                
            Lapl.apply( src, dst, level, flag );
            }
        }  

        void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                  mat,
                  const typename VecFunctionType::template FunctionType< idx_t >&   src,
                  const typename VecFunctionType::template FunctionType< idx_t >&   dst,
                  size_t                                                            level,
                  DoFType                                                           flag ) const
        {
            if(nu > 0) {
                WALBERLA_ABORT("AL currently unsupported");
     
         
            // Sparse matrix objects
            PETScSparseMatrix<LaplOpType> LaplMat;
            PETScSparseMatrix<DivOpType> divMat;
            PETScSparseMatrix<DivOpTType> divTMat;
         
            // numerators
            typename LaplOpType::srcType::template FunctionType< idx_t > VelocityNum( "VelocityNum", storage_, level_, level_ );
            typename DivOpType::dstType::template FunctionType< idx_t > PressureNum( "PressureNum", storage_, level_, level_ );
            VelocityNum.enumerate(level_); 
            PressureNum.enumerate(level_);
            
            // assembly
            LaplMat.createMatrixFromOperator(Lapl, level_, VelocityNum);
            divMat.createMatrixFromOperator(div, level_, VelocityNum, PressureNum);
            divTMat.createMatrixFromOperator(divT, level_, PressureNum, VelocityNum);
            auto divTdivMat = mat->createCopy();
            
            // Proxy objects for createFromLinearCombination
            auto LaplMatProxy = std::make_shared<PETScSparseMatrixProxy>(LaplMat.get());
            auto divMatProxy = std::make_shared<PETScSparseMatrixProxy>(divMat.get());
            auto divTMatProxy = std::make_shared<PETScSparseMatrixProxy>(divTMat.get());

            // final Augmented Lagrangian matrix creation
            divTdivMat->createFromMatrixProduct({divTMatProxy,divMatProxy});
            mat->createFromMatrixLinComb({1,nu},{LaplMatProxy,divTdivMat});
            
            } else {
                Lapl.toMatrix(mat,src,dst,level,flag);
            }
        }   
    private:
        std::shared_ptr< PrimitiveStorage >          storage_;
        size_t                                       level_;
        real_t                                       nu;
        LaplOpType                                   Lapl;
        DivOpType                                    div;
        DivOpTType                                   divT;
};

// specification of Augmented Lagrangian Operator for a P2P1 Taylor Hood Stokes discretization
using ALOP_P2P1TH = AugmentedLagrangianOperator< 
    P2VectorFunction<real_t>, 
    P2ConstantVectorLaplaceOperator, 
    P2ToP1ConstantDivOperator, 
    P1ToP2ConstantDivTOperator 
>;

} // namespace hyteg
