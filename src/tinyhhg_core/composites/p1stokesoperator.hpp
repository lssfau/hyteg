#pragma once

#include "tinyhhg_core/composites/p1stokesfunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"

#include "tinyhhg_core/p1functionspace/P1HelperFunctions.hpp"

namespace hhg
{

class P1StokesOperator
{
public:

  P1StokesOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel, std::shared_ptr<std::array<P1Function<real_t>*, 2>> normals = {})
    : A(storage, minLevel, maxLevel),
      div_x(storage, minLevel, maxLevel),
      div_y(storage, minLevel, maxLevel),
      divT_x(storage, minLevel, maxLevel),
      divT_y(storage, minLevel, maxLevel),
      pspg(storage, minLevel, maxLevel),
      normals_(normals),
      storage_(storage)

  {
    freeslip_ = normals != nullptr;
  }

  void apply(P1StokesFunction<real_t>& src, P1StokesFunction<real_t>& dst, size_t level, DoFType flag)
  {
    DoFType newFlag = flag;

    if (isFreeslip()) {
      newFlag = newFlag | DirichletBoundary;
    }

    A.apply(src.u, dst.u, level, newFlag, Replace);
    divT_x.apply(src.p, dst.u, level, newFlag, Add);

    A.apply(src.v, dst.v, level, newFlag, Replace);
    divT_y.apply(src.p, dst.v, level, newFlag, Add);

    if (isFreeslip()) {
      // project normals
      P1::projectNormal<2>(storage_, {{&src.u, &src.v}}, *normals_, level, DirichletBoundary);
    }

    div_x.apply(src.u, dst.p, level, newFlag | DirichletBoundary, Replace);
    div_y.apply(src.v, dst.p, level, newFlag | DirichletBoundary, Add);
    pspg.apply(src.p, dst.p, level, newFlag | DirichletBoundary, Add);
  }

  P1LaplaceOperator A;
  P1DivxOperator div_x;
  P1DivyOperator div_y;
  P1DivTxOperator divT_x;
  P1DivTyOperator divT_y;
  P1PSPGOperator pspg;

  bool isFreeslip() const {
    return freeslip_;
  }

  const std::shared_ptr<std::array<P1Function<real_t>*, 2>>& getNormals() const {
    return normals_;
  };

private:
  bool freeslip_;
  std::shared_ptr<std::array<P1Function<real_t>*, 2>> normals_;
  const std::shared_ptr< PrimitiveStorage > & storage_;
};

}