
template<>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Primitive >( const PrimitiveID & id ) const { return primitiveExistsLocally( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Vertex >( const PrimitiveID & id ) const { return vertexExistsLocally( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Edge >  ( const PrimitiveID & id ) const { return edgeExistsLocally( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Face >  ( const PrimitiveID & id ) const { return faceExistsLocally( id ); }


template<>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Primitive >( const PrimitiveID & id ) const { return primitiveExistsInNeighborhood( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Vertex >( const PrimitiveID & id ) const { return vertexExistsInNeighborhood( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Edge >  ( const PrimitiveID & id ) const { return edgeExistsInNeighborhood( id ); }

template<>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Face >  ( const PrimitiveID & id ) const { return faceExistsInNeighborhood( id ); }


template<>
inline const Primitive* PrimitiveStorage::getPrimitiveGenerically< Primitive >( const PrimitiveID & id ) const { return getPrimitive( id ); }

template<>
inline       Primitive* PrimitiveStorage::getPrimitiveGenerically< Primitive >( const PrimitiveID & id )       { return getPrimitive( id ); }

template<>
inline const Vertex* PrimitiveStorage::getPrimitiveGenerically< Vertex >( const PrimitiveID & id ) const { return getVertex( id ); }

template<>
inline       Vertex* PrimitiveStorage::getPrimitiveGenerically< Vertex >( const PrimitiveID & id )       { return getVertex( id ); }

template<>
inline const Edge*   PrimitiveStorage::getPrimitiveGenerically< Edge >  ( const PrimitiveID & id ) const { return getEdge( id ); }

template<>
inline       Edge*   PrimitiveStorage::getPrimitiveGenerically< Edge >  ( const PrimitiveID & id )       { return getEdge( id ); }

template<>
inline const Face*   PrimitiveStorage::getPrimitiveGenerically< Face >  ( const PrimitiveID & id ) const { return getFace( id ); }

template<>
inline       Face*   PrimitiveStorage::getPrimitiveGenerically< Face >  ( const PrimitiveID & id )       { return getFace( id ); }

template<>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Vertex >( std::vector< PrimitiveID > & primitiveIDs ) const { getVertexIDs( primitiveIDs ); }

template<>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Edge >( std::vector< PrimitiveID > & primitiveIDs ) const { getEdgeIDs( primitiveIDs ); }

template<>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Face >( std::vector< PrimitiveID > & primitiveIDs ) const { getFaceIDs( primitiveIDs ); }
