#ifndef __ML_MESH_TYPE_H
#define __ML_MESH_TYPE_H

#include <vector>
#include <algorithm/complex_base.h>

// Forward declarations needed for creating the used types
class CVertexO;
class CEdgeO;
class CFaceO;

// Declaration of the semantic of the used types
class CUsedTypesO: public tessellation::UsedTypes < tessellation::Use<CVertexO>::AsVertexType,
    tessellation::Use<CFaceO  >::AsFaceType >{};


// The Main Vertex Class
// Most of the attributes are optional and must be enabled before use.
// Each vertex needs 40 byte, on 32bit arch. and 44 byte on 64bit arch.

class CVertexO  : public tessellation::Vertex< CUsedTypesO,
    tessellation::vertex::InfoOcf,           /*  4b */
    tessellation::vertex::Coord3f,           /* 12b */
    tessellation::vertex::BitFlags,          /*  4b */
    tessellation::vertex::Normal3f,          /* 12b */
    tessellation::vertex::VFAdjOcf,          /*  0b */
    tessellation::vertex::TexCoordfOcf      /*  0b */
>{    
};

// Each face needs 32 byte, on 32bit arch. and 48 byte on 64bit arch.
class CFaceO : public tessellation::Face<  CUsedTypesO,
    tessellation::face::InfoOcf,              /* 4b */
    tessellation::face::VertexRef,  /*12b */
    tessellation::face::BitFlags,             /* 4b */
    tessellation::face::Normal3f,             /*12b */
    tessellation::face::FFAdjOcf,             /* 0b */
    tessellation::face::VFAdjOcf,             /* 0b */
    tessellation::face::MarkOcf
> {
};

class CMeshO  : public tessellation::tri::TriMesh<tessellation::vertex::vector_ocf<CVertexO>, tessellation::face::vector_ocf<CFaceO> >
{
};

#endif
