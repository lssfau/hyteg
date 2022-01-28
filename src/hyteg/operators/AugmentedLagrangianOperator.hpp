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
#include "hyteg/operators/Operator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/mixedoperators/P1ScalarToP2VectorOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include <type_traits>
namespace hyteg {

using walberla::real_t;

template < class VecFunctionType, class LaplOpType, class DivOpType, class DivOpTType > // different discretizations for Laplace operator and Divergence
class AugmentedLagrangianOperator : public Operator<  VecFunctionType,  VecFunctionType >
{

    public:
    
        AugmentedLagrangianOperator(
            const std::shared_ptr< PrimitiveStorage >& storage,
            size_t                                     level, // only one level, GKB not to be leveled any time soon
            real_t                                     gamma // gamma = Augmented Laplacian Parameter
        ) 
        : Operator< VecFunctionType, VecFunctionType >( storage, level, level )
        , storage_(storage)
        , level_(level)
        , Lapl( storage, level, level )
        , div( storage, level, level )
        , divT( storage, level, level )
        { 
            gamma_ = gamma;
        }

        void apply( 
            const VecFunctionType&         src,
            const VecFunctionType&         dst,
            size_t                     level,
            DoFType                    flag,
            UpdateType                 updateType = Replace 
        ) const
        {
            // apply Augmented Lagrangian opterator:
            // ALOPx = (Laplace + gamma^-1 Div^T Div)x
            // in [AR13] notation: M = W + gamma^-1AA^T
            VecFunctionType tmp0("tmp0",storage_, level, level);
            typename DivOpType::dstType tmp1("tmp1",storage_, level, level);
            VecFunctionType tmp2("tmp2",storage_, level, level);
            Lapl.apply( src, tmp0, level, flag );
            if(fabs(gamma_) > 1e-10) {
                div.apply( src, tmp1, level, flag );
                divT.apply( tmp1, tmp2, level, flag );
                dst.assign( {1, 1/gamma_}, {tmp0, tmp2}, level );
            } else {
                // if gamma is 0 switch off AL
                dst.assign( {1}, {tmp0}, level );   
            }
        }  

        void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                  mat,
                  const typename VecFunctionType::template FunctionType< idx_t >&   src,
                  const typename VecFunctionType::template FunctionType< idx_t >&   dst,
                  size_t                                                            level,
                  DoFType                                                           flag ) const
   {
        /*
        auto divMat = mat->createCopy();
        auto divTMat = mat->createCopy();
        
        div.toMatrix(divMat,src,dst,level,flag);
        divT.toMatrix(divTMat,src,dst,level,flag);
        mat->createFromMatrixProduct({divTMat,divMat});
        */

       // use gamma = 0 for now
       //TOTO create AL matrix
       Lapl.toMatrix(mat,src,dst,level,flag);

   } 
    private:
        std::shared_ptr< PrimitiveStorage >          storage_;
        size_t                                       level_;
        real_t                                       gamma_;
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
