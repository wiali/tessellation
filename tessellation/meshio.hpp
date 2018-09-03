#ifndef MESHLOADING_HPP_
#define MESHLOADING_HPP_

#include "ml_mesh_type.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <tuple>
#include <QVector3D>
#include <QVector2D>

using namespace std;

namespace MeshSpace
{
    inline std::vector<std::string> split(
        const std::string& str,
        const char delimiter
    ) {
        std::vector<std::string> internal;
        std::stringstream ss(str);
        std::string tok;

        while (getline(ss, tok, delimiter)) {
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
                return -1;

            mesh.vert.EnableTexCoord();

            size_t nVertexs = raw_vertices.size();
            CMeshO::VertexIterator vi = tessellation::tri::Allocator<CMeshO>::AddVertices(mesh, nVertexs);

            size_t nFaces = indices.size() / 3;
            CMeshO::FaceIterator fi = tessellation::tri::Allocator<CMeshO>::AddFaces(mesh, nFaces);


            CMeshO::VertexPointer *ivp = new CMeshO::VertexPointer[nVertexs];
            for (int i=0; i < nVertexs; ++i)
            {
                QVector3D pt = raw_vertices[i];
                ivp[i] = &*vi; vi->P() = CMeshO::CoordType(pt[0], pt[1], pt[2]); ++vi;
            }

            CMeshO::FacePointer *ifp = new CMeshO::FacePointer[nFaces];
            for (int i = 0; i < nFaces; ++i)
            {
                int nIndex1 = indices[i * 3];
                int nIndex2 = indices[i * 3 + 1];
                int nIndex3 = indices[i * 3 + 2];
                ifp[i] = &*fi; 
                fi->V(0) = ivp[nIndex1];
                fi->V(1) = ivp[nIndex2];
                fi->V(2) = ivp[nIndex3];

                fi->V(0)->T().u() = uvs[nIndex1][0];
                fi->V(0)->T().v() = uvs[nIndex1][1];
                fi->V(0)->T().n() = 0;

                fi->V(1)->T().u() = uvs[nIndex2][0];
                fi->V(1)->T().v() = uvs[nIndex2][1];
                fi->V(1)->T().n() = 0;

                fi->V(2)->T().u() = uvs[nIndex3][0];
                fi->V(2)->T().v() = uvs[nIndex3][1];
                fi->V(2)->T().n() = 0;

                fi->V(0)->N().Import(tessellation::Point3f(normals[nIndex1][0], normals[nIndex1][1], normals[nIndex1][2]));
                fi->V(1)->N().Import(tessellation::Point3f(normals[nIndex2][0], normals[nIndex2][1], normals[nIndex2][2]));
                fi->V(2)->N().Import(tessellation::Point3f(normals[nIndex3][0], normals[nIndex3][1], normals[nIndex3][2]));

                ++fi;

            }


            return 0;
        }

        static int writeMesh(std::string path, CMeshO& mesh)
        {
            std::ofstream fs;
            fs.open(path, std::ofstream::out);

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

#endif
