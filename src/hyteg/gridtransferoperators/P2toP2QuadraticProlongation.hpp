#pragma once

#include "core/DataTypes.h"

#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"

namespace hyteg {

template < typename VType >
class P2Function;

class P2toP2QuadraticProlongation : public ProlongationOperator< P2Function< walberla::real_t > >
{
 public:
   void prolongate( const P2Function< walberla::real_t >& function,
                    const walberla::uint_t&               sourceLevel,
                    const DoFType&                        flag ) const override;

 private:
   void prolongateAdditively( const P2Function< walberla::real_t >& function,
                              const walberla::uint_t&               sourceLevel,
                              const DoFType&                        flag ) const;

   void prolongateAdditively3D( const P2Function< walberla::real_t >& function,
                                const walberla::uint_t&               sourceLevel,
                                const DoFType&                        flag ) const;

   void prolongateStandard( const P2Function< walberla::real_t >& function,
                            const walberla::uint_t&               sourceLevel,
                            const DoFType&                        flag ) const;
};

} // namespace hyteg
