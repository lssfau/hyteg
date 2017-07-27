#include "p1bubblefunction.hpp"
#include "tinyhhg_core/Function.hpp"
#include "p1bubblevertex.hpp"
#include "p1bubbleedge.hpp"
#include "p1bubbleface.hpp"

namespace hhg {

P1BubbleFunction::P1BubbleFunction(const std::string& name,
                                   const std::shared_ptr< PrimitiveStorage > & storage,
                                   uint_t minLevel,
                                   uint_t maxLevel)
    : Function(name, storage, minLevel, maxLevel)
{
  // TODO: add P1BubbleDataHandling (to be implemented)
}

P1BubbleFunction::~P1BubbleFunction()
{
    //TODO implement!
}

void P1BubbleFunction::interpolate(std::function<real_t(const hhg::Point3D&)>& expr, uint_t level, DoFType flag = All)
{

}

void P1BubbleFunction::assign(const std::vector<walberla::real_t> scalars, const std::vector<P1BubbleFunction*> functions, size_t level, DoFType flag = All)
{

}

void P1BubbleFunction::add(const std::vector<walberla::real_t> scalars, const std::vector<P1BubbleFunction*> functions, size_t level, DoFType flag = All)
{

}

real_t P1BubbleFunction::dot(P1BubbleFunction& rhs, size_t level, DoFType flag = All)
{

}

void P1BubbleFunction::prolongate(size_t level, DoFType flag = All){

}

void P1BubbleFunction::prolongateQuadratic(size_t level, DoFType flag = All){

}

void P1BubbleFunction::restrict(size_t level, DoFType flag = All)
{

}
}
