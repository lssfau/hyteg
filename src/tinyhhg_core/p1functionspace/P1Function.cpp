#include "tinyhhg_core/p1functionspace/P1Function.hpp"

namespace hhg
{

//FIXME remove after we are in walberla namespace
using namespace walberla::mpistubs;


  P1Function::P1Function(const std::string& name, const PrimitiveStorage& storage, uint_t minLevel, uint_t maxLevel) :
  Function(name,storage,minLevel,maxLevel)
  {
    FaceP1FunctionMemoryDataHandling faceP1FunctionMemoryDataHandling;
    EdgeP1FunctionMemoryDataHandling edgeP1FunctionMemoryDataHandling;
    VertexP1FunctionMemoryDataHandling vertexP1FunctionMemoryDataHandling;
    faceDataID_   = storage->addFaceData(faceP1FunctionMemoryDataHandling,name);
    edgeDataID_   = storage->addEdgeData(edgeP1FunctionMemoryDataHandling,name);
    vertexDataID_ = storage->addVertexData(vertexP1FunctionMemoryDataHandling,name);
  }

  ~P1Function();

  void interpolate(std::function<real_t(const hhg::Point3D&)>& expr, uint_t level, DoFType flag = All);

  void assign(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag = All);

  void add(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag = All);

  real_t dot(P1Function& rhs, size_t level, DoFType flag = All);

  void prolongate(size_t level, DoFType flag = All);

  void prolongateQuadratic(size_t level, DoFType flag = All);

  void restrict(size_t level, DoFType flag = All);

}
