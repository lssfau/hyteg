
#pragma once

namespace hhg {

class P1Transport
{
public:

  P1Transport( const std::shared_ptr< PrimitiveStorage > & storage, const uint_t minLevel, const uint_t maxLevel ) :
    invLumpedMass_( storage, minLevel, maxLevel ), A_( storage, minLevel, maxLevel ),
    divT_x_( storage, minLevel, maxLevel ), divT_y_( storage, minLevel, maxLevel ), divT_z_( storage, minLevel, maxLevel ),
    tmp0_( "tmp0", storage, minLevel, maxLevel ),
    tmp1_( "tmp1", storage, minLevel, maxLevel ),
    tmp2_( "tmp2", storage, minLevel, maxLevel ),
    tmp3_( "tmp3", storage, minLevel, maxLevel ),
    tmp4_( "tmp4", storage, minLevel, maxLevel )
  {}

  void step( const P1Function< real_t > & c, 
             const P1Function< real_t > & ux, const P1Function< real_t > & uy, const P1Function< real_t > & uz,
             const uint_t & level, const DoFType & flag, const real_t & dt, const real_t & a )
  {
    WALBERLA_ASSERT_GREATER_EQUAL( level, minLevel );
    WALBERLA_ASSERT_GREATER_EQUAL( maxLevel, level );

    A_.apply( c, tmp0_, level, flag, Replace );
    tmp0_.assign( {a * dt}, {&tmp0_}, level, flag );

    divT_x_.apply( c, tmp1_, level, flag, Replace );
    divT_y_.apply( c, tmp2_, level, flag, Replace );
    divT_z_.apply( c, tmp3_, level, flag, Replace );
    tmp1_.multElementwise( {&ux, &tmp1_}, level, flag );
    tmp2_.multElementwise( {&uy, &tmp2_}, level, flag );
    tmp3_.multElementwise( {&uz, &tmp3_}, level, flag );

    tmp4_.assign( {dt, dt, dt}, {&tmp1_, &tmp2_, &tmp3_}, level, flag );

    tmp0_.add( {1.0}, {&tmp4_}, level, flag );
    invLumpedMass_.apply( tmp0_, c, level, flag, Add );
  }
  

private:

  P1InvLumpedMassOperator   invLumpedMass_;
  P1ConstantLaplaceOperator A_;
  P1DivTxOperator           divT_x_;
  P1DivTyOperator           divT_y_;
  P1DivTzOperator           divT_z_;
  P1Function< real_t >      tmp0_;
  P1Function< real_t >      tmp1_;
  P1Function< real_t >      tmp2_;
  P1Function< real_t >      tmp3_;
  P1Function< real_t >      tmp4_;

}

}


