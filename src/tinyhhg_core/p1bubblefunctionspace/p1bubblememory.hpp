#pragma once

#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include "tinyhhg_core/levelinfo.hpp"

#include <string>


namespace hhg
{

  //******************************************* MEMORY CLASSES **************************************//

  class VertexP1BubbleStencilMemory
    :public VertexMemory
  {
  public:
    VertexP1BubbleStencilMemory() : VertexMemory(Stencil) { ; }

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

    ~VertexP1BubbleStencilMemory() { free(); }

  };


  class VertexP1BubbleMemory
    :public VertexMemory
  {
  public:
    VertexP1BubbleMemory() : VertexMemory(P1) { ; }

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

    ~VertexP1BubbleMemory() { free(); }
  };


  class EdgeP1BubbleStencilMemory
    :public EdgeMemory
  {
  public:
    EdgeP1BubbleStencilMemory() : EdgeMemory(Stencil) { ; }

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
      //WALBERLA_LOG_DEVEL("EdgeStencilMemory, kind = " + std::to_string(this->type));
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

    ~EdgeP1BubbleStencilMemory() { free(); }

  };


  class EdgeP1BubbleMemory
    :public EdgeMemory
  {
  public:
    EdgeP1BubbleMemory() : EdgeMemory(P1) { ; }

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

    ~EdgeP1BubbleMemory() { free(); }
  };


  class FaceP1BubbleStencilMemory
    :public FaceMemory
  {
  public:
    FaceP1BubbleStencilMemory() : FaceMemory(Stencil) { ; }

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

    ~FaceP1BubbleStencilMemory() { free(); }

  };


  class FaceP1BubbleMemory
    :public FaceMemory
  {
  public:
    FaceP1BubbleMemory() : FaceMemory(P1) { ; }

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
      return levelinfo::num_microfaces_per_face(level);
    }


    ~FaceP1BubbleMemory() { free(); }
  };



  //******************************************* MEMORY CASTS ****************************************//


  inline VertexP1BubbleMemory* getVertexP1BubbleMemory(Vertex& vertex, size_t id)
  {
#ifndef NDEBUG
    if (vertex.memory.size() <= id)
      WALBERLA_LOG_WARNING("Memory ID is out of range");
      if (vertex.memory[id]->type != P1)
        WALBERLA_LOGLEVEL_WARNING("Trying to convert something to VertexP1Memory which is not of the right type");
#endif // !

    return static_cast<VertexP1BubbleMemory*>(vertex.memory[id]);
  }

  inline VertexP1BubbleStencilMemory* getVertexP1BubbleStencilMemory(const Vertex& vertex, size_t id)
  {
    return static_cast<VertexP1BubbleStencilMemory*>(vertex.memory[id]);
  }

  inline EdgeP1BubbleStencilMemory* getEdgeP1BubbleStencilMemory(const Edge& edge, size_t id)
  {
    return static_cast<EdgeP1BubbleStencilMemory*>(edge.memory[id]);
  }

  inline EdgeP1BubbleMemory* getEdgeP1BubbleMemory(const Edge& edge, size_t id)
  {
    return static_cast<EdgeP1BubbleMemory*>(edge.memory[id]);
  }

  inline FaceP1BubbleStencilMemory* getFaceP1BubbleStencilMemory(const Face& face, size_t id)
  {
    return static_cast<FaceP1BubbleStencilMemory*>(face.memory[id]);
  }

  inline FaceP1BubbleMemory* getFaceP1BubbleMemory(const Face& face, size_t id)
  {
    return static_cast<FaceP1BubbleMemory*>(face.memory[id]);
  }
}