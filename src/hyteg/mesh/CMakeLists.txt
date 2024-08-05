target_sources( hyteg
    PRIVATE
    GmshReaderForMSH22.hpp
    GmshReaderForMSH41.cpp
    GmshReaderForMSH41.hpp
    HyTeGMeshDir.cpp
    MeshGenAnnulus.cpp
    MeshGenCubedDomain.cpp
    MeshGenCuboid.cpp
    MeshGenFaceChain.cpp
    MeshGenFromGmshFile.cpp
    MeshGenRectangle.cpp
    MeshGenRefineCoarseMesh.cpp
    MeshGenSingleSimplex.cpp
    MeshGenSphericalShell.cpp
    MeshGenSymmetricCuboid.cpp     
    MeshGenTorus.cpp
    MeshGenUnitSquare.cpp
    MeshInfo.cpp
    MeshInfo.hpp
)

add_subdirectory(micro)