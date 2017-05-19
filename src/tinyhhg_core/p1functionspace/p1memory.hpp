#ifndef P1MEMORY_HPP
#define P1MEMORY_HPP

#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include "tinyhhg_core/levelinfo.hpp"


namespace hhg
{

  //******************************************* MEMORY CLASSES **************************************//

  class VertexStencilMemory
    :public VertexMemory
  {
  public:
    VertexStencilMemory() : VertexMemory(Stencil) { ; }

    std::map<size_t, double*> data;
    size_t num_deps;

    inline virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline double* addlevel(size_t level, size_t num_deps)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        this->num_deps = num_deps;
        data[level] = new double[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      return levelinfo::num_microvertices_per_vertex(level) + num_deps;
    }

    ~VertexStencilMemory() { free(); }

  };


  class VertexP1Memory
    :public VertexMemory
  {
  public:
    VertexP1Memory() : VertexMemory(P1) { ; }

    std::map<size_t, double*> data;
    size_t num_deps;

    inline virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline double* addlevel(size_t level, size_t num_deps)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        this->num_deps = num_deps;
        data[level] = new double[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      return levelinfo::num_microvertices_per_vertex(level) + num_deps;
    }

    ~VertexP1Memory() { free(); }
  };


  class EdgeStencilMemory
    :public EdgeMemory
  {
  public:
    EdgeStencilMemory() : EdgeMemory(Stencil) { ; }

    std::map<size_t, double*> data;

    virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline double* addlevel(size_t level)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        data[level] = new double[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      return 7;
    }

    ~EdgeStencilMemory() { free(); }

  };


  class EdgeP1Memory
    :public EdgeMemory
  {
  public:
    EdgeP1Memory() : EdgeMemory(P1) { ; }

    std::map<size_t, double*> data;
    size_t num_deps;

    virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline double* addlevel(size_t level, size_t num_deps)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        this->num_deps = num_deps;
        data[level] = new double[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      size_t num_dofs_per_edge = levelinfo::num_microvertices_per_edge(level);
      return num_dofs_per_edge + num_deps*(num_dofs_per_edge-1);
    }

    ~EdgeP1Memory() { free(); }
  };


  class FaceStencilMemory
    :public FaceMemory
  {
  public:
    FaceStencilMemory() : FaceMemory(Stencil) { ; }

    std::map<size_t, double*> data;

    virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline double* addlevel(size_t level)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        data[level] = new double[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      return 7;
    }

    ~FaceStencilMemory() { free(); }

  };


  class FaceP1Memory
    :public FaceMemory
  {
  public:
    FaceP1Memory() : FaceMemory(P1) { ; }

    std::map<size_t, double*> data;

    virtual void free()
    {
      for (auto el : data)
      {
        delete[] el.second;
      }
      data.clear();
    }

    inline double* addlevel(size_t level)
    {
      if (data.count(level)>0)
        WALBERLA_LOG_WARNING("Level already exists.")
      else
      {
        data[level] = new double[getSize(level)]();
      }
      return data[level];
    }

    inline size_t getSize(size_t level)
    {
      return levelinfo::num_microfaces_per_face(level);
    }


    ~FaceP1Memory() { free(); }
  };



  //******************************************* MEMORY CASTS ****************************************//


  inline VertexP1Memory* getVertexP1Memory(Vertex& vertex, size_t id)
  {
#ifndef NDEBUG
    if (vertex.memory.size() <= id)
      WALBERLA_LOG_WARNING("Memory ID is out of range");
      if (vertex.memory[id]->type != P1)
        WALBERLA_LOGLEVEL_WARNING("Trying to convert something to VertexP1Memory which is not of the right type");
#endif // !

    return static_cast<VertexP1Memory*>(vertex.memory[id]);
  }

  inline VertexStencilMemory* getVertexStencilMemory(const Vertex& vertex, size_t id)
  {
    return static_cast<VertexStencilMemory*>(vertex.memory[id]);
  }

  inline EdgeStencilMemory* getEdgeStencilMemory(const Edge& edge, size_t id)
  {
    return static_cast<EdgeStencilMemory*>(edge.memory[id]);
  }

  inline EdgeP1Memory* getEdgeP1Memory(const Edge& edge, size_t id)
  {
    return static_cast<EdgeP1Memory*>(edge.memory[id]);
  }

  inline FaceStencilMemory* getFaceStencilMemory(const Face& face, size_t id)
  {
    return static_cast<FaceStencilMemory*>(face.memory[id]);
  }

  inline FaceP1Memory* getFaceP1Memory(const Face& face, size_t id)
  {
    return static_cast<FaceP1Memory*>(face.memory[id]);
  }
}


#endif // !P1MEMORY_HPP
