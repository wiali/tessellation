/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/


#ifndef __VCG_TRIMESH_CLOSEST
#define __VCG_TRIMESH_CLOSEST
#include <math.h>

#include <algorithm/space_base.h>
#include <algorithm/box3.h>
#include <algorithm/grid_closest.h>
#include <algorithm/simplex_base.h>
#include <algorithm/complex_base.h>

namespace tessellation {
    namespace tri {

        //**MARKER CLASSES**//
        template <class MESH_TYPE, class OBJ_TYPE>
        class Tmark
        {
            MESH_TYPE *m;
        public:
            Tmark() {}
            Tmark(MESH_TYPE *m) { SetMesh(m); }
            void UnMarkAll() { tessellation::tri::UnMarkAll(*m); }
            bool IsMarked(OBJ_TYPE* obj) { return (tessellation::tri::IsMarked(*m, obj)); }
            void Mark(OBJ_TYPE* obj) { tessellation::tri::Mark(*m, obj); }
            void SetMesh(MESH_TYPE *_m)
            {
                tri::RequirePerFaceMark(*_m);
                m = _m;
            }
        };

        template <class MESH_TYPE>
        class FaceTmark :public Tmark<MESH_TYPE, typename MESH_TYPE::FaceType>
        {
        public:
            FaceTmark() {}
            FaceTmark(MESH_TYPE *m) { this->SetMesh(m); }
        };

        template <class MESH_TYPE>
        class EdgeTmark :public Tmark<MESH_TYPE, typename MESH_TYPE::EdgeType>
        {
        public:
            EdgeTmark() {}
            EdgeTmark(MESH_TYPE *m) { this->SetMesh(m); }
        };


        template <class MESH_TYPE>
        class EmptyTMark
        {
        public:
            typedef typename  MESH_TYPE::VertexType VertexType;
            typedef typename  MESH_TYPE::EdgeType EdgeType;
            typedef typename  MESH_TYPE::FaceType FaceType;
            inline EmptyTMark() {}
            inline EmptyTMark(MESH_TYPE *) {}
            inline void UnMarkAll() const {}
            inline bool IsMarked(VertexType*) const { return false; }
            inline void Mark(VertexType*) const {}
            inline bool IsMarked(EdgeType*) const { return false; }
            inline void Mark(EdgeType*) const {}
            inline bool IsMarked(FaceType*) const { return false; }
            inline void Mark(FaceType*) const {}
            inline void SetMesh(void * /*m=0*/) const {}
        };

        //**CLOSEST FUNCTION DEFINITION**//

        /*

        aka MetroCore
        data una mesh m e una ug sulle sue facce trova il punto di m piu' vicino ad
        un punto dato.
        */

        // input: mesh, punto, griglia (gr), distanza limite (mdist)
        // output: normale (interpolata) alla faccia e punto piu' vicino su di essa, e coord baricentriche del punto trovato

        // Nota che il parametro template GRID non ci dovrebbe essere, visto che deve essere
        // UGrid<MESH::FaceContainer >, ma non sono riuscito a definirlo implicitamente

        template <class MESH, class GRID>
        typename MESH::FaceType * GetClosestFaceEP(MESH & mesh, GRID & gr, const typename GRID::CoordType & _p,
            const typename GRID::ScalarType & _maxDist, typename GRID::ScalarType & _minDist,
            typename GRID::CoordType & _closestPt, typename GRID::CoordType & _normf,
            typename GRID::CoordType & _ip)
        {
            typedef typename GRID::ScalarType ScalarType;

            typedef FaceTmark<MESH> MarkerFace;
            MarkerFace mf(&mesh);
            vcg::face::PointDistanceEPFunctor<ScalarType> FDistFunct;
            _minDist = _maxDist;
            typename MESH::FaceType* bestf = gr.GetClosest(FDistFunct, mf, _p, _maxDist, _minDist, _closestPt);

            if (_maxDist > ScalarType(fabs(_minDist)))
            {
                // f=bestf;
                //calcolo normale con interpolazione trilineare
                InterpolationParameters<typename MESH::FaceType, typename MESH::ScalarType>(*bestf, bestf->N(), _closestPt, _ip);
                _normf = (bestf->V(0)->cN())*_ip[0] +
                    (bestf->V(1)->cN())*_ip[1] +
                    (bestf->V(2)->cN())*_ip[2];

                _minDist = fabs(_minDist);
                return(bestf);
            }
            return (0);
        }

        template <class MESH, class GRID>
        typename MESH::FaceType * GetClosestFaceBase(MESH & mesh, GRID & gr, const typename GRID::CoordType & _p,
            const typename GRID::ScalarType _maxDist, typename GRID::ScalarType & _minDist,
            typename GRID::CoordType &_closestPt)
        {
            typedef typename GRID::ScalarType ScalarType;
            typedef FaceTmark<MESH> MarkerFace;
            MarkerFace mf;
            mf.SetMesh(&mesh);
            vcg::face::PointDistanceBaseFunctor<ScalarType> PDistFunct;
            _minDist = _maxDist;
            return (gr.GetClosest(PDistFunct, mf, _p, _maxDist, _minDist, _closestPt));
        }

        template <class MESH, class GRID>
        typename MESH::FaceType * GetClosestFaceBase(MESH & mesh, GRID & gr, const typename GRID::CoordType & _p,
            const typename GRID::ScalarType _maxDist, typename GRID::ScalarType & _minDist,
            typename GRID::CoordType &_closestPt, typename GRID::CoordType & _normf,
            typename GRID::CoordType & _ip)
        {
            typedef typename GRID::ScalarType ScalarType;
            typename MESH::FaceType * f = GetClosestFaceBase(mesh, gr, _p, _maxDist, _minDist, _closestPt);
            if (_maxDist > ScalarType(fabs(_minDist)))
            {
                // normal computed with trilinear interpolation
                InterpolationParameters<typename MESH::FaceType, typename MESH::ScalarType>(*f, f->N(), _closestPt, _ip);
                _normf = (f->V(0)->cN())*_ip[0] +
                    (f->V(1)->cN())*_ip[1] +
                    (f->V(2)->cN())*_ip[2];
            }
            return f;
        }

        template <class MESH, class GRID>
        typename MESH::FaceType * GetClosestFaceEP(MESH & mesh, GRID & gr, const typename GRID::CoordType & _p,
            const typename GRID::ScalarType _maxDist, typename GRID::ScalarType & _minDist,
            typename GRID::CoordType &_closestPt)
        {
            typedef typename GRID::ScalarType ScalarType;
            typedef FaceTmark<MESH> MarkerFace;
            MarkerFace mf;
            mf.SetMesh(&mesh);
            vcg::face::PointDistanceEPFunctor<ScalarType> PDistFunct;
            _minDist = _maxDist;
            return (gr.GetClosest(PDistFunct, mf, _p, _maxDist, _minDist, _closestPt));
        }

        template <class MESH, class GRID>
        typename MESH::FaceType * GetClosestFaceNormal(MESH & mesh, GRID & gr, const typename MESH::VertexType & _p,
            const typename GRID::ScalarType & _maxDist, typename GRID::ScalarType & _minDist,
            typename GRID::CoordType &_closestPt)
        {
            typedef FaceTmark<MESH> MarkerFace;
            MarkerFace mf;
            mf.SetMesh(&mesh);
            typedef vcg::face::PointNormalDistanceFunctor<typename MESH::VertexType> PDistFunct;
            PDistFunct fn;
            _minDist = _maxDist;
            //return (gr.GetClosest(PDistFunct,mf,_p,_maxDist,_minDist,_closestPt.P()));
            return (gr.template GetClosest <PDistFunct, MarkerFace>(fn, mf, _p, _maxDist, _minDist, _closestPt));
        }

        template <class MESH, class GRID>
        typename MESH::VertexType * GetClosestVertex(MESH & mesh, GRID & gr, const typename GRID::CoordType & _p,
            const typename GRID::ScalarType & _maxDist, typename GRID::ScalarType & _minDist)
        {
            typedef typename GRID::ScalarType ScalarType;
            typedef Point3<ScalarType> Point3x;
            typedef EmptyTMark<MESH> MarkerVert;
            MarkerVert mv;
            mv.SetMesh(&mesh);
            typedef vcg::vertex::PointDistanceFunctor<typename MESH::ScalarType> VDistFunct;
            VDistFunct fn;
            _minDist = _maxDist;
            Point3x _closestPt;
            return (gr.template GetClosest<VDistFunct, MarkerVert>(fn, mv, _p, _maxDist, _minDist, _closestPt));
        }

        template <class MESH, class GRID>
        typename MESH::VertexType * GetClosestVertexScale(MESH & mesh, GRID & gr, const typename GRID::CoordType & _p,
            const typename MESH::CoordType & center, const typename GRID::ScalarType & _maxDist, typename GRID::ScalarType & _minDist)
        {
            typedef typename GRID::ScalarType ScalarType;
            typedef Point3<ScalarType> Point3x;
            typedef EmptyTMark<MESH> MarkerVert;
            MarkerVert mv;
            mv.SetMesh(&mesh);
            typedef vcg::vertex::PointScaledDistanceFunctor<typename MESH::VertexType> VDistFunct;
            VDistFunct fn;
            fn.Cen() = center;
            _minDist = _maxDist;
            Point3x _closestPt;
            return (gr.template GetClosest<VDistFunct, MarkerVert>(fn, mv, _p, _maxDist, _minDist, _closestPt));
        }

        template <class MESH, class GRID>
        typename MESH::VertexType * GetClosestVertexNormal(MESH & mesh, GRID & gr, const typename MESH::VertexType & _p,
            const typename GRID::ScalarType & _maxDist, typename GRID::ScalarType & _minDist)
        {
            typedef typename GRID::ScalarType ScalarType;
            typedef Point3<ScalarType> Point3x;
            typedef EmptyTMark<MESH> MarkerVert;
            MarkerVert mv;
            mv.SetMesh(&mesh);
            typedef vcg::vertex::PointNormalDistanceFunctor<typename MESH::VertexType> VDistFunct;
            VDistFunct fn;
            _minDist = _maxDist;
            Point3x _closestPt;
            return (gr.template GetClosest <VDistFunct, MarkerVert>(fn, mv, _p, _maxDist, _minDist, _closestPt));
        }

        template <class MESH, class GRID, class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
        unsigned int GetKClosestFaceEP(MESH & mesh, GRID & gr, const unsigned int _k,
            const typename GRID::CoordType & _p, const typename GRID::ScalarType & _maxDist,
            OBJPTRCONTAINER & _objectPtrs, DISTCONTAINER & _distances, POINTCONTAINER & _points)
        {
            typedef FaceTmark<MESH> MarkerFace;
            MarkerFace mf;
            mf.SetMesh(&mesh);
            vcg::face::PointDistanceEPFunctor<typename MESH::ScalarType> FDistFunct;
            return (gr.GetKClosest /*<vcg::face::PointDistanceFunctor, MarkerFace,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>*/
            (FDistFunct, mf, _k, _p, _maxDist, _objectPtrs, _distances, _points));
        }

        template <class MESH, class GRID, class OBJPTRCONTAINER>
        unsigned int GetInBoxFace(MESH & mesh,
            GRID & gr,
            const tessellation::Box3<typename GRID::ScalarType> _bbox,
            OBJPTRCONTAINER & _objectPtrs)
        {
            typedef FaceTmark<MESH> MarkerFace;
            MarkerFace mf(&mesh);
            return(gr.GetInBox/*<MarkerFace,OBJPTRCONTAINER>*/(mf, _bbox, _objectPtrs));
        }

        template <class MESH, class GRID, class OBJPTRCONTAINER>
        unsigned int GetInBoxVertex(MESH & mesh,
            GRID & gr,
            const tessellation::Box3<typename GRID::ScalarType> _bbox,
            OBJPTRCONTAINER & _objectPtrs)
        {
            typedef EmptyTMark<MESH> MarkerVert;
            MarkerVert mv;
            mv.SetMesh(&mesh);
            return(gr.GetInBox/*<MarkerVert,OBJPTRCONTAINER>*/(mv, _bbox, _objectPtrs));
        }

    }	 // end namespace tri
}	 // end namespace vcg

#endif
