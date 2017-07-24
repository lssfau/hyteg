#pragma once
#include "tinyhhg_core/Function.hpp"

namespace hhg
{

class VertexP1FunctionMemory;
class EdgeP1FunctionMemory;
class FaceP1FunctionMemory;

class P1Function : public Function
{
public:

  P1Function(const std::string& name, PrimitiveStorage& storage, uint_t minLevel, uint_t maxLevel);

  ~P1Function();

  void interpolate(std::function<real_t(const hhg::Point3D&)>& expr, uint_t level, DoFType flag = All);

  void assign(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag = All);

  void add(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag = All);

  real_t dot(P1Function& rhs, size_t level, DoFType flag = All);

  void prolongate(size_t level, DoFType flag = All);

  void prolongateQuadratic(size_t level, DoFType flag = All);

  void restrict(size_t level, DoFType flag = All);

private:
  PrimitiveDataID<VertexP1FunctionMemory,Vertex> vertexDataID_;
  PrimitiveDataID<EdgeP1FunctionMemory  ,Edge  > edgeDataID_;
  PrimitiveDataID<FaceP1FunctionMemory  ,Face  > faceDataID_;

};

}
