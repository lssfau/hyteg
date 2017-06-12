#ifndef P1MEMORY_HPP
#define P1MEMORY_HPP

#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include "tinyhhg_core/levelinfo.hpp"

#include <string>


namespace hhg
{

class VertexP1StencilMemory
  :public VertexMemory
{
public:
  VertexP1StencilMemory() : VertexMemory(MemoryType::P1Stencil) { ; }

  std::map<size_t, real_t*> data;
  size_t num_deps;

  inline virtual void free()
  {
    for (auto el : data)
    {
      delete[] el.second;
    }
    data.clear();
  }

  inline real_t* addlevel(size_t level, size_t num_deps)
  {
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      this->num_deps = num_deps;
      data[level] = new real_t[getSize(level)]();
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    return levelinfo::num_microvertices_per_vertex(level) + num_deps;
  }

  ~VertexP1StencilMemory() { free(); }

};


class VertexP1FunctionMemory
  :public VertexMemory
{
public:
  VertexP1FunctionMemory() : VertexMemory(MemoryType::P1Function) { ; }

  std::map<size_t, real_t*> data;
  size_t num_deps;

  inline virtual void free()
  {
    for (auto el : data)
    {
      delete[] el.second;
    }
    data.clear();
  }

  inline real_t* addlevel(size_t level, size_t num_deps)
  {
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      this->num_deps = num_deps;
      data[level] = new real_t[getSize(level)]();
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    return levelinfo::num_microvertices_per_vertex(level) + num_deps;
  }

  ~VertexP1FunctionMemory() { free(); }
};


class EdgeP1StencilMemory
  :public EdgeMemory
{
public:
  EdgeP1StencilMemory() : EdgeMemory(MemoryType::P1Stencil) { ; }

  std::map<size_t, real_t*> data;

  virtual void free()
  {
    for (auto el : data)
    {
      delete[] el.second;
    }
    data.clear();
  }

  inline real_t* addlevel(size_t level)
  {
    //WALBERLA_LOG_DEVEL("EdgeP1StencilMemory, kind = " + std::to_string(this->type));
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      data[level] = new real_t[getSize(level)]();
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    return 7;
  }

  ~EdgeP1StencilMemory() { free(); }

};


class EdgeP1FunctionMemory
  :public EdgeMemory
{
public:
  EdgeP1FunctionMemory() : EdgeMemory(MemoryType::P1Function) { ; }

  std::map<size_t, real_t*> data;
  size_t num_deps;

  virtual void free()
  {
    for (auto el : data)
    {
      delete[] el.second;
    }
    data.clear();
  }

  inline real_t* addlevel(size_t level, size_t num_deps)
  {
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      this->num_deps = num_deps;
      data[level] = new real_t[getSize(level)]();
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    size_t num_dofs_per_edge = levelinfo::num_microvertices_per_edge(level);
    return num_dofs_per_edge + num_deps*(num_dofs_per_edge-1);
  }

  ~EdgeP1FunctionMemory() { free(); }
};


class FaceP1StencilMemory
  :public FaceMemory
{
public:
  FaceP1StencilMemory() : FaceMemory(MemoryType::P1Stencil) { ; }

  std::map<size_t, real_t*> data;

  virtual void free()
  {
    for (auto el : data)
    {
      delete[] el.second;
    }
    data.clear();
  }

  inline real_t* addlevel(size_t level)
  {
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      data[level] = new real_t[getSize(level)]();
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    return 7;
  }

  ~FaceP1StencilMemory() { free(); }

};


class FaceP1FunctionMemory
  :public FaceMemory
{
public:
  FaceP1FunctionMemory() : FaceMemory(MemoryType::P1Function) { ; }

  std::map<size_t, real_t*> data;

  virtual void free()
  {
    for (auto el : data)
    {
      delete[] el.second;
    }
    data.clear();
  }

  inline real_t* addlevel(size_t level)
  {
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      data[level] = new real_t[getSize(level)]();
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    return levelinfo::num_microvertices_per_face(level);
  }


  ~FaceP1FunctionMemory() { free(); }
};

namespace P1 {

inline VertexP1StencilMemory *getVertexStencilMemory(const Vertex &vertex, size_t id) {
  return static_cast<VertexP1StencilMemory *>(vertex.memory[id]);
}

inline VertexP1FunctionMemory *getVertexFunctionMemory(Vertex &vertex, size_t id) {
#ifndef NDEBUG
  if (vertex.memory.size() <= id) WALBERLA_LOG_WARNING("Memory ID is out of range");
  if (vertex.memory[id]->type != MemoryType::P1Function)
      WALBERLA_LOGLEVEL_WARNING("Trying to convert something to VertexP1FunctionMemory which is not of the right type");
#endif // !

  return static_cast<VertexP1FunctionMemory *>(vertex.memory[id]);
}

inline EdgeP1StencilMemory *getEdgeStencilMemory(const Edge &edge, size_t id) {
  return static_cast<EdgeP1StencilMemory *>(edge.memory[id]);
}

inline EdgeP1FunctionMemory *getEdgeFunctionMemory(const Edge &edge, size_t id) {
  return static_cast<EdgeP1FunctionMemory *>(edge.memory[id]);
}

inline FaceP1StencilMemory *getFaceStencilMemory(const Face &face, size_t id) {
  return static_cast<FaceP1StencilMemory *>(face.memory[id]);
}

inline FaceP1FunctionMemory *getFaceFunctionMemory(const Face &face, size_t id) {
  return static_cast<FaceP1FunctionMemory *>(face.memory[id]);
}

}
}


#endif // !P1MEMORY_HPP
