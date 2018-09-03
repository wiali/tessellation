#include "stdafx.h"
#include "meshio.hpp"
#include <algorithm/complex_base.h>
#include <algorithm/split_triangle.h>
#include <QImage>

using namespace tessellation;

enum MeshElement
{
    MM_NONE = 0x00000000,
    MM_VERTFACETOPO = 0x00000040,
    MM_VERTTEXCOORD = 0x00000400,
    MM_FACEMARK = 0x00020000,
    MM_FACEFACETOPO = 0x00040000,
    MM_ALL = 0xffffffff
};

void updateDataMask(CMeshO& mesh, int neededDataMask)
{
    if ((neededDataMask & MM_FACEFACETOPO) != 0)
    {
        mesh.face.EnableFFAdjacency();
        tri::UpdateTopology<CMeshO>::FaceFace(mesh);
    }
    if ((neededDataMask & MM_VERTFACETOPO) != 0)
    {
        mesh.vert.EnableVFAdjacency();
        mesh.face.EnableVFAdjacency();
        tri::UpdateTopology<CMeshO>::VertexFace(mesh);
    }

    if ((neededDataMask & MM_FACEMARK) != 0)
        mesh.face.EnableMark();
    if ((neededDataMask & MM_VERTTEXCOORD) != 0)
        mesh.vert.EnableTexCoord();
}

void UpdateBoxAndNormals(CMeshO& mesh)
{
    tri::UpdateNormal<CMeshO>::PerFaceNormalized(mesh);
    tri::UpdateNormal<CMeshO>::PerVertexAngleWeighted(mesh);
    for (auto i = mesh.vert.begin(); i != mesh.vert.end(); i++)
    {
        i->N() = i->N().Normalize();
    }
}


class GetClosestFace
{
    typedef tessellation::GridStaticPtr<CMeshO::FaceType, CMeshO::ScalarType > MetroMeshGrid;
    typedef tri::FaceTmark<CMeshO> MarkerFace;

public:

    GetClosestFace() {}

    void init(CMeshO *_m)
    {
        m = _m;
        if (m)
        {
            unifGrid.Set(m->face.begin(), m->face.end());
            markerFunctor.SetMesh(m);
            dist_upper_bound = 10.0;// m->bbox.Diag() / 10.0f;
        }
    }

    CMeshO *m;

    MetroMeshGrid unifGrid;

    MarkerFace markerFunctor;

    float dist_upper_bound;

    CMeshO::FaceType * getFace(Point3f &p)
    {
        assert(m);
        // the results
        Point3f closestPt;
        float dist = dist_upper_bound;
        const CMeshO::CoordType &startPt = p;

        // compute distance between startPt and the mesh S2
        CMeshO::FaceType   *nearestF = 0;
        tessellation::face::PointDistanceBaseFunctor<CMeshO::ScalarType> PDistFunct;
        dist = dist_upper_bound;

        nearestF = unifGrid.GetClosest(PDistFunct, markerFunctor, startPt, dist_upper_bound, dist, closestPt);

        //if (dist == dist_upper_bound) qDebug() << "Dist is = upper bound";

        return nearestF;
    }
};

/** evaluate barycentric coordinates
@param bq Point on the face
@param L0 barycentric value for V(0)
@param L1 barycentric value for V(1)
@param L2 barycentric value for V(2)
@return true se bq appartain to the face, false otherwise
from http://en.wikipedia.org/wiki/Barycentric_coordinate_system_(mathematics)
L1=((y2-y3)(x-x3)+(x3-x2)(y-y3))/((y2-y3)(x1-x3)+(x3-x2)(y1-y3))
L2=((y3-y1)(x-x3)+(x1-x3)(y-y3))/((y3-y1)(x2-x3)+(x1-x3)(y2-y3))
L3=1-L1-L2
*/
Point3f Barycentric(Point3f aV1, Point3f aV2, Point3f aV3, Point3f aP)
{
    Point3f a = aV2 - aV3, b = aV1 - aV3, c = aP - aV3;
    float aLen = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    float bLen = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    float ab = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    float ac = a[0] * c[0] + a[1] * c[1] + a[2] * c[2];
    float bc = b[0] * c[0] + b[1] * c[1] + b[2] * c[2];
    float d = aLen * bLen - ab * ab;

    Point3f ret;
    ret[0] = (aLen * bc - ab * ac) / d;
    ret[1] = (bLen * ac - ab * bc) / d;
    ret[2] = 1.0f - ret[0] - ret[1];

    return ret;
}

void updateUVs(CMeshO& mesh, CMeshO& mesh_old)
{
    auto pGetClosestFace = new GetClosestFace();
    updateDataMask(mesh_old, MM_FACEMARK);
    //set up the 
    pGetClosestFace->init(&mesh_old);

    for (auto i = mesh.vert.begin(); i != mesh.vert.end(); i++)
    {
        float u = i->T().P()[0];
        float v = i->T().P()[1];

        //https://answers.unity.com/questions/1374972/calculate-uv-coordinates-based-on-known-vertices-a.html
        //if (math::isZero(u) && math::isZero(v))
        if (fabs(u) > 1 && fabs(v) > 1)
        {
            auto pFace = pGetClosestFace->getFace(i->P());
            Point3f fBarycentric;
            fBarycentric = Barycentric(pFace->V(0)->P(), pFace->V(1)->P(), pFace->V(2)->P(), i->P());
            auto newUV = pFace->V(0)->T().P() * fBarycentric[0] + pFace->V(1)->T().P() * fBarycentric[1] + pFace->V(2)->T().P() * fBarycentric[2];
            i->T().P()[0] = newUV[0];
            i->T().P()[1] = newUV[1];
        }
    }
}


static int copyMesh(CMeshO& mesh_src, CMeshO& mesh_dest)
{
    //mesh_dest.vert.EnableTexCoord();

    mesh_dest.vert.EnableTexCoord();

    size_t nVertexs = mesh_src.vert.size();
    CMeshO::VertexIterator vi = tessellation::tri::Allocator<CMeshO>::AddVertices(mesh_dest, nVertexs);

    size_t nFaces = mesh_src.face.size();
    CMeshO::FaceIterator fi = tessellation::tri::Allocator<CMeshO>::AddFaces(mesh_dest, nFaces);

    CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[nVertexs];
    for (int i = 0; i < nVertexs; ++i)
    {
        ivp[i] = &*vi;
        vi->P() = mesh_src.vert[i].P();
        vi->T().P()[0] = mesh_src.vert[i].T().P()[0];
        vi->T().P()[1] = mesh_src.vert[i].T().P()[1];
        vi->N() = mesh_src.vert[i].N();
        ++vi;
    }

    CMeshO::FacePointer *ifp = new CMeshO::FacePointer[nFaces];
    for (int i = 0; i < nFaces; ++i)
    {
        size_t nIndex1 = tessellation::tri::Index(mesh_src, mesh_src.face[i].V(0));
        size_t nIndex2 = tessellation::tri::Index(mesh_src, mesh_src.face[i].V(1));
        size_t nIndex3 = tessellation::tri::Index(mesh_src, mesh_src.face[i].V(2));
        ifp[i] = &*fi;
        fi->V(0) = ivp[nIndex1];
        fi->V(1) = ivp[nIndex2];
        fi->V(2) = ivp[nIndex3];

        fi->V(0)->T().u() = mesh_src.vert[nIndex1].T().P()[0];
        fi->V(0)->T().v() = mesh_src.vert[nIndex1].T().P()[1];
        fi->V(0)->T().n() = 0;

        fi->V(1)->T().u() = mesh_src.vert[nIndex2].T().P()[0];
        fi->V(1)->T().v() = mesh_src.vert[nIndex2].T().P()[1];
        fi->V(1)->T().n() = 0;

        fi->V(2)->T().u() = mesh_src.vert[nIndex3].T().P()[0];
        fi->V(2)->T().v() = mesh_src.vert[nIndex3].T().P()[1];
        fi->V(2)->T().n() = 0;

        fi->V(0)->N().Import(mesh_src.vert[nIndex1].N());
        fi->V(1)->N().Import(mesh_src.vert[nIndex2].N());
        fi->V(2)->N().Import(mesh_src.vert[nIndex3].N());

        fi->N() = tri::TriangleNormal(*fi).Normalize();

        ++fi;

    }

    tri::Allocator<CMeshO>::CompactFaceVector(mesh_dest);
    tri::Allocator<CMeshO>::CompactVertexVector(mesh_dest);
    updateDataMask(mesh_dest, MM_FACEFACETOPO);
    tri::UpdateFlags<CMeshO>::FaceBorderFromFF(mesh_dest);

    updateDataMask(mesh_dest, MM_VERTTEXCOORD);
    updateDataMask(mesh_dest, MM_VERTFACETOPO);

    return 0;
}

int main()
{
    CMeshO mesh;

    CMeshO mesh_old;

    MeshSpace::MeshIO<MeshSpace::MeshFileType::OBJ>::loadMesh("agate.obj", mesh);

    copyMesh(mesh, mesh_old);

    tri::Allocator<CMeshO>::CompactFaceVector(mesh);
    tri::Allocator<CMeshO>::CompactVertexVector(mesh);
    updateDataMask(mesh, MM_FACEFACETOPO);
    tri::UpdateFlags<CMeshO>::FaceBorderFromFF(mesh);

    updateDataMask(mesh, MM_VERTTEXCOORD);

    int iterations = 1;
    for (int i = 0; i < iterations; ++i)
    {
        updateDataMask(mesh, MM_VERTFACETOPO);

        RefineOddEven(mesh, tri::OddPointLoop<CMeshO>(mesh), tri::EvenPointLoop<CMeshO>(), 0.0);
    }

    UpdateBoxAndNormals(mesh);

    updateUVs(mesh, mesh_old);

    QImage imageHeight;
    imageHeight.load("Height.jpg");
        
    int nHeightMapSize = 2048;
    float heightMapScale = 2.0;
    for (auto i = mesh.vert.begin(); i != mesh.vert.end(); i++)
    {
        auto pt = i->P();
        auto uv = i->T().P();
        auto _normal = i->N();

        int _x = uv[0] * nHeightMapSize;
        int _y = uv[1] * nHeightMapSize;
        
        float pixelValue = qRed(imageHeight.pixel(_x, _y)) / 255.0;

        float displacement = pixelValue;
        i->P() += _normal * displacement * heightMapScale;
    }

    UpdateBoxAndNormals(mesh);

    updateUVs(mesh, mesh_old);

    MeshSpace::MeshIO<MeshSpace::MeshFileType::OBJ>::writeMesh("A_dog.obj", mesh);

    return 0;
}

