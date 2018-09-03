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

///** Returns the closest posistion of a point _p and its distance
//@param _p a 3d point
//@param _maxDist maximum distance not to search beyond.
//@param _getPointDistance (templated type) a functor object used to calculate distance from a grid object to the point _p.
//@param _minDist the returned closest distance
//@param _closestPt the returned closest point
//@return The closest element
//*/
///*
//	A DISTFUNCT object must implement an operator () with signature:
//		bool operator () (const ObjType& obj, const CoordType & _p, ScalarType & min_dist, CoordType & _closestPt);
//*/

#pragma once

#include <algorithm/space_base.h>

namespace tessellation {

    template <class SPATIAL_INDEX, class OBJPOINTDISTFUNCTOR, class OBJMARKER>
    typename SPATIAL_INDEX::ObjPtr  GridClosest(SPATIAL_INDEX &Si,
        OBJPOINTDISTFUNCTOR _getPointDistance,
        OBJMARKER & _marker,
        const typename OBJPOINTDISTFUNCTOR::QueryType  & _p_obj,
        const typename SPATIAL_INDEX::ScalarType & _maxDist,
        typename SPATIAL_INDEX::ScalarType & _minDist,
        typename SPATIAL_INDEX::CoordType &_closestPt)
    {
        typedef typename SPATIAL_INDEX::ObjPtr ObjPtr;
        typedef typename SPATIAL_INDEX::CoordType CoordType;
        typedef typename SPATIAL_INDEX::ScalarType ScalarType;
        typedef typename SPATIAL_INDEX::Box3x Box3x;

        Point3<ScalarType> _p = OBJPOINTDISTFUNCTOR::Pos(_p_obj);
        // Initialize min_dist with _maxDist to exploit early rejection test.
        _minDist = _maxDist;

        ObjPtr winner = NULL;
        _marker.UnMarkAll();
        ScalarType newradius = Si.voxel.Norm();
        ScalarType radius;
        Box3i iboxdone, iboxtodo;
        CoordType t_res;
        typename SPATIAL_INDEX::CellIterator first, last, l;
        if (Si.bbox.IsInEx(_p))
        {
            Point3i _ip;
            Si.PToIP(_p, _ip);
            Si.Grid(_ip[0], _ip[1], _ip[2], first, last);
            for (l = first; l != last; ++l)
            {
                ObjPtr elem = &(**l);
                if (!elem->IsD())
                {
                    if (_getPointDistance((**l), _p_obj, _minDist, t_res))  // <-- NEW: use of distance functor
                    {
                        winner = elem;
                        _closestPt = t_res;
                        newradius = _minDist; //
                    }
                    _marker.Mark(elem);
                }
            }
            iboxdone = Box3i(_ip, _ip);
        }

        int ix, iy, iz;
        Box3i ibox(Point3i(0, 0, 0), Si.siz - Point3i(1, 1, 1));
        do
        {
            radius = newradius;
            Box3x boxtodo = Box3x(_p, radius);
            //boxtodo.Intersect(Si.bbox);
            Si.BoxToIBox(boxtodo, iboxtodo);
            iboxtodo.Intersect(ibox);
            if (!boxtodo.IsNull())
            {
                for (ix = iboxtodo.min[0]; ix <= iboxtodo.max[0]; ix++)
                    for (iy = iboxtodo.min[1]; iy <= iboxtodo.max[1]; iy++)
                        for (iz = iboxtodo.min[2]; iz <= iboxtodo.max[2]; iz++)
                            if (ix<iboxdone.min[0] || ix>iboxdone.max[0] ||  // this test is to avoid to re-process already analyzed cells.
                                iy<iboxdone.min[1] || iy>iboxdone.max[1] ||
                                iz<iboxdone.min[2] || iz>iboxdone.max[2])
                            {
                                Si.Grid(ix, iy, iz, first, last);
                                for (l = first; l != last; ++l) if (!(**l).IsD())
                                {
                                    ObjPtr elem = &(**l);
                                    if (!elem->IsD())
                                    {
                                        if (!_marker.IsMarked(elem))
                                        {
                                            if (_getPointDistance((**l), _p_obj, _minDist, t_res))
                                            {
                                                winner = elem;
                                                _closestPt = t_res;
                                            };
                                            _marker.Mark(elem);
                                        }
                                    }
                                }
                            }
            }
            if (!winner) newradius = radius + Si.voxel.Norm();
            else newradius = _minDist;
            iboxdone = iboxtodo;
        } while (_minDist > radius);

        return winner;
    }

    template <class SPATIALINDEXING, class OBJPOINTDISTFUNCTOR, class OBJMARKER,
        class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
        unsigned int GridGetKClosest(SPATIALINDEXING &_Si,
            OBJPOINTDISTFUNCTOR & _getPointDistance,
            OBJMARKER & _marker,
            const unsigned int _k,
            const typename SPATIALINDEXING::CoordType & _p,
            const typename SPATIALINDEXING::ScalarType & _maxDist,
            OBJPTRCONTAINER & _objectPtrs,
            DISTCONTAINER & _distances,
            POINTCONTAINER & _points)
    {
        typedef vcg::ClosestIterator<SPATIALINDEXING, OBJPOINTDISTFUNCTOR, OBJMARKER> ClosestIteratorType;
        ClosestIteratorType	Cli = ClosestIteratorType(_Si, _getPointDistance);
        Cli.SetMarker(_marker);
        Cli.Init(_p, _maxDist);
        unsigned int i = 0;
        _objectPtrs.clear();
        _distances.clear();
        _points.clear();
        while ((!Cli.End()) && (i < _k))
        {
            _objectPtrs.push_back(&(*Cli));
            _distances.push_back(Cli.Dist());
            _points.push_back(Cli.NearestPoint());
            ++Cli;
            i++;
        }
        return (i);
    }

}//end namespace vcg

