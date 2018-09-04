#pragma once

#include <complex_base.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <tuple>
#include <QVector3D>
#include <QVector2D>
#include <QFileInfo>

using namespace std;
using namespace tessellation;

namespace MeshSpace
{
    inline std::vector<std::string> split(
        const std::string& str,
        const char delimiter
    )
    {
        std::vector<std::string> internal;
        std::stringstream ss(str);
        std::string tok;
        while (getline(ss, tok, delimiter))
        {
            internal.push_back(tok);
        }
        return internal;
    }

    enum MeshFileType
    {
        OBJ, OFF
    };

    template<MeshFileType FileType>
    struct MeshIO
    {

        static int loadMesh(
            std::string,
            std::vector<QVector3D>&,
            std::vector<int>&
        )
        {
            std::cout << "Not implemented.\n";
            return -1;
        }

        static int loadMesh(std::string)
        {
            std::cout << "Not implemented.\n";
            return -1;
        }
    };

    template<>
    struct MeshIO<MeshFileType::OBJ>
    {
        static int loadMesh(
            std::string path,
            std::vector<QVector3D>& vertices,
            std::vector<int>& indices,
            std::vector<QVector2D>& uvs,
            std::vector<QVector3D>& normals
        )
        {
            std::ifstream fs;
            fs.open(path, std::ifstream::in);
            if (!fs)
            {
                std::cout << "Mesh file \"" << path << "\" could not be loaded." << std::endl;
                return -1;
            }
            std::vector<QVector3D> single_normals;
            std::vector<int> normal_indices;
            std::vector<QVector2D> single_uvs;
            std::vector<int> uvs_indices;
            while (!fs.eof())
            {
                std::string type;
                fs >> type;
                if (type[0] == '#')
                {
                    std::getline(fs, type);;
                }
                if (type == "v")
                {
                    QVector3D new_vertex;
                    fs >> new_vertex[0];
                    fs >> new_vertex[1];
                    fs >> new_vertex[2];
                    vertices.push_back(new_vertex);
                }
                if (type == "vn")
                {
                    QVector3D new_normal;
                    fs >> new_normal[0];
                    fs >> new_normal[1];
                    fs >> new_normal[2];
                    single_normals.push_back(new_normal);
                }
                if (type == "vt")
                {
                    QVector2D new_uv;
                    fs >> new_uv[0];
                    fs >> new_uv[1];
                    single_uvs.push_back(new_uv);
                }
                if (type == "f")
                {
                    std::string vertex_desc;
                    for (int i = 0; i < 3; ++i)
                    {
                        fs >> vertex_desc;
                        std::vector<std::string> vert_split;
                        vert_split = split(vertex_desc, '/');
                        indices.push_back(std::stoi(vert_split[0]) - 1);
                    }
                }
            }
            uvs = single_uvs;
            normals = single_normals;
            return 0;
        }

        static int loadMesh(std::string path, CMeshO& mesh)
        {
            std::vector<QVector3D> raw_vertices;
            std::vector<int> indices;
            std::vector<QVector2D> uvs;
            std::vector<QVector3D> normals;
            if (loadMesh(path, raw_vertices, indices, uvs, normals) < 0)
            {
                return -1;
            }
            
            return loadMesh(raw_vertices, indices, uvs, normals, mesh);
        }

        // For vinnie MD structure data convert
        static int loadMesh(const std::vector<float>& vertices, const std::vector<int>& indices,
            const std::vector<float>& uvs, const std::vector<float>& normals, const QMatrix4x4& transform, CMeshO& mesh)
        {
            std::vector<QVector3D> _vertices;
            std::vector<int> _indices;
            std::vector<QVector2D> _uvs;
            std::vector<QVector3D> _normals;

            _vertices.reserve(vertices.size());
            _indices.reserve(indices.size());
            _uvs.reserve(uvs.size());
            _normals.reserve(normals.size());

            for (int i = 0; i < vertices.size(); i += 3)
                _vertices.push_back(transform * QVector3D(vertices[i], vertices[i + 1], vertices[i + 2]));

            for (int i = 0; i < uvs.size(); i += 2)
                _uvs.push_back(QVector2D(uvs[i], 1.0f - uvs[i + 1]));

            for (int i = 0; i < normals.size(); i += 3)
            {
                auto pt = transform * QVector3D(normals[i], normals[i + 1], normals[i + 2]);
                pt.normalize();
                _normals.push_back(pt);
            }

            for (int i = 0; i < indices.size(); i += 9)
            {
                _indices.push_back(indices[i]);
                _indices.push_back(indices[i + 3]);
                _indices.push_back(indices[i + 6]);
            }

            loadMesh(_vertices, _indices, _uvs, _normals, mesh);
        }

        static int loadMesh(const std::vector<QVector3D>& raw_vertices, const std::vector<int>& indices, 
            const std::vector<QVector2D>& uvs, const std::vector<QVector3D>& normals, CMeshO& mesh)
        {
            mesh.vert.EnableTexCoord();
            size_t nVertexs = raw_vertices.size();
            CMeshO::VertexIterator vi = tessellation::tri::Allocator<CMeshO>::AddVertices(mesh, nVertexs);
            size_t nFaces = indices.size() / 3;
            CMeshO::FaceIterator fi = tessellation::tri::Allocator<CMeshO>::AddFaces(mesh, nFaces);
            CMeshO::VertexPointer* ivp = new CMeshO::VertexPointer[nVertexs];
            for (int i = 0; i < nVertexs; ++i)
            {
                QVector3D pt = raw_vertices[i];
                ivp[i] = &*vi;
                vi->P() = CMeshO::CoordType(pt[0], pt[1], pt[2]);
                vi->T().u() = uvs[i][0];
                vi->T().v() = uvs[i][1];
                vi->T().n() = 0;
                vi->N() = tessellation::Point3f(normals[i][0], normals[i][1], normals[i][2]);
                ++vi;
            }
            CMeshO::FacePointer* ifp = new CMeshO::FacePointer[nFaces];
            for (int i = 0; i < nFaces; ++i)
            {
                int nIndex1 = indices[i * 3];
                int nIndex2 = indices[i * 3 + 1];
                int nIndex3 = indices[i * 3 + 2];
                ifp[i] = &*fi;
                fi->V(0) = ivp[nIndex1];
                fi->V(1) = ivp[nIndex2];
                fi->V(2) = ivp[nIndex3];
                ++fi;
            }
            return 0;
        }

        static int writeMesh(std::string path, CMeshO& mesh)
        {
            std::ofstream fs;
            fs.open(path + ".obj", std::ofstream::out);

            string basename = QFileInfo(path.c_str()).baseName().toStdString();

            fs << "mtllib " << basename << ".mtl" << endl;
            fs << "usemtl material_0" << endl;
            fs << endl;

            for (auto i = mesh.vert.begin(); i != mesh.vert.end(); i++)
            {
                fs << "v " << i->P().X() << " " << i->P().Y() << " " << i->P().Z() << std::endl;
                auto _normal = i->N();
                fs << "vn " << _normal[0] << " " << _normal[1] << " " << _normal[2] << std::endl;
                fs << "vt " << i->T().P()[0] << " " << i->T().P()[1] << std::endl;
            }
            size_t arrIndex[3];
            for (auto i = mesh.face.begin(); i != mesh.face.end(); i++)
            {
                fs << "f ";
                // +1 because Obj file format begins from index = 1 but not from index = 0.
                arrIndex[0] = tessellation::tri::Index(mesh, (*i).V(0)) + 1;
                arrIndex[1] = tessellation::tri::Index(mesh, (*i).V(1)) + 1;
                arrIndex[2] = tessellation::tri::Index(mesh, (*i).V(2)) + 1;
                fs << arrIndex[0] << "/" << arrIndex[0] << "/" << arrIndex[0] << " "
                   << arrIndex[1] << "/" << arrIndex[1] << "/" << arrIndex[1] << " "
                   << arrIndex[2] << "/" << arrIndex[2] << "/" << arrIndex[2] << endl;
            }
            return 0;
        }
    };
}
