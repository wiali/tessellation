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

#pragma once

#include <vector>
#include <algorithm>
#include <stdio.h>


namespace tessellation
{

    /** BasicGrid

    Basic Class abstracting a gridded structure in a 3d space;
    Usueful for having coherent float to integer conversion in a unique place:
    Some Notes:
    - bbox is the real occupation of the box in the space;
    - siz is the number of cells for each side

    OBJTYPE:      Type of the indexed objects.
    SCALARTYPE:   Scalars type for structure's internal data (may differ from
    object's scalar type).

    */

    template <class SCALARTYPE>
    class BasicGrid
    {
    public:

        typedef SCALARTYPE ScalarType;
        typedef Box3<ScalarType> Box3x;
        typedef Point3<ScalarType> CoordType;
        typedef BasicGrid<SCALARTYPE> GridType;

        Box3x bbox;

        CoordType dim;		/// Spatial Dimention (edge legth) of the bounding box
        Point3i siz;		/// Number of cells forming the grid
        CoordType voxel;	/// Dimensions of a single cell


        /* Given a 3D point p, returns the index of the corresponding cell
        * @param p is a 3D point in the space
        * @return integer coordinates pi of the cell
        */
        inline void PToIP(const CoordType & p, Point3i &pi) const
        {
            CoordType t = p - bbox.min;
            pi[0] = int(t[0] / voxel[0]);
            pi[1] = int(t[1] / voxel[1]);
            pi[2] = int(t[2] / voxel[2]);
        }

        /* Given a cell in <ScalarType> coordinates, compute the corresponding cell in integer coordinates
        * @param b is the cell in <ScalarType> coordinates
        * @return ib is the correspondent box in integer coordinates
        */
        inline void BoxToIBox(const Box3x & b, Box3i & ib) const
        {
            PToIP(b.min, ib.min);
            PToIP(b.max, ib.max);
            //assert(ib.max[0]>=0 && ib.max[1]>=0 && ib.max[2]>=0);
        }
    };


    /** Calcolo dimensioni griglia.
    Calcola la dimensione della griglia in funzione
    della ratio del bounding box e del numero di elementi
    */
    template<class scalar_type>
    void BestDim(const __int64 elems, const Point3<scalar_type> & size, Point3i & dim)
    {
        const __int64 mincells = 1;		// Numero minimo di celle
        const double GFactor = 1;	// GridEntry = NumElem*GFactor
        double diag = size.Norm();	// Diagonale del box
        double eps = diag*1e-4;		// Fattore di tolleranza

        assert(elems>0);
        assert(size[0] >= 0.0);
        assert(size[1] >= 0.0);
        assert(size[2] >= 0.0);


        __int64 ncell = (__int64)(elems*GFactor);	// Calcolo numero di voxel
        if (ncell<mincells)
            ncell = mincells;

        dim[0] = 1;
        dim[1] = 1;
        dim[2] = 1;

        if (size[0]>eps)
        {
            if (size[1]>eps)
            {
                if (size[2]>eps)
                {
                    double k = pow((double)(ncell / (size[0] * size[1] * size[2])), double(1.0 / 3.f));
                    dim[0] = int(size[0] * k);
                    dim[1] = int(size[1] * k);
                    dim[2] = int(size[2] * k);
                }
                else
                {
                    dim[0] = int(::sqrt(ncell*size[0] / size[1]));
                    dim[1] = int(::sqrt(ncell*size[1] / size[0]));
                }
            }
            else
            {
                if (size[2]>eps)
                {
                    dim[0] = int(::sqrt(ncell*size[0] / size[2]));
                    dim[2] = int(::sqrt(ncell*size[2] / size[0]));
                }
                else
                    dim[0] = int(ncell);
            }
        }
        else
        {
            if (size[1]>eps)
            {
                if (size[2]>eps)
                {
                    dim[1] = int(::sqrt(ncell*size[1] / size[2]));
                    dim[2] = int(::sqrt(ncell*size[2] / size[1]));
                }
                else
                    dim[1] = int(ncell);
            }
            else if (size[2]>eps)
                dim[2] = int(ncell);
        }
        dim[0] = std::max(dim[0], 1);
        dim[1] = std::max(dim[1], 1);
        dim[2] = std::max(dim[2], 1);
    }

    /****************************************************************************
    Class SpatialIndex

    Description:
    This class exposes the base interface for all spatial indexing data
    structures, i.e. grids, bounding volume trees.

    Template Parameters:
    OBJTYPE:      Type of the indexed objects.
    SCALARTYPE:   Scalars type for structure's internal data (may differ from
    object's scalar type).

    ****************************************************************************/

    template <class OBJTYPE, class SCALARTYPE>
    class SpatialIndex {
    public:
        /**************************************************************************
        Commonly used typedefs.
        **************************************************************************/
        typedef SpatialIndex<OBJTYPE, SCALARTYPE> ClassType;
        typedef OBJTYPE ObjType;
        typedef SCALARTYPE ScalarType;
        typedef ObjType * ObjPtr;
        typedef Point3<ScalarType> CoordType;
        typedef Box3<ScalarType> BoxType;

        
        /**************************************************************************
        Method GetInBox.

        Description:
        The GetInBox returns all the object in the specified bbox

        Template Parameters:

        OBJMARKER           : The type of a marker functor.
        OBJPTRCONTAINER     : The type of a object pointers container.

        Method Parameters:
        _marker           : [IN] Functor for marking objects already tested.
        _bbox             : [IN] The bounding box of spatial query.
        _objectPtrs       : [OUT] Container which, in return, will contain pointers to the closest objects.


        Return Value:
        The number of in-box objects found.

        **************************************************************************/
        template <class OBJMARKER, class OBJPTRCONTAINER>
        unsigned int GetInBox(OBJMARKER & _marker, const BoxType _bbox, OBJPTRCONTAINER & _objectPtrs) {
            assert(0);
            (void)_marker;
            (void)_bbox;
            (void)_objectPtrs;
            return (0);
        }

    };

    /** Static Uniform Grid
    A spatial search structure for a accessing a container of objects.
    It is based on a uniform grid overlayed over a protion of space.
    The grid partion the space into cells. Cells contains just pointers
    to the object that are stored elsewhere.
    The set of objects is meant to be static and pointer stable.

    Useful for situation were many space related query are issued over
    the same dataset (ray tracing, measuring distances between meshes,
    re-detailing ecc.).
    Works well for distribution that ar reasonably uniform.
    How to use it:
    ContainerType must have a 'value_type' typedef inside.
    (stl containers already have it)

    Objects pointed by cells (of kind 'value_type') must have
    a 'ScalarType' typedef (float or double usually)
    and a member function:

    void GetBBox(Box3<ScalarType> &b)
    which return the bounding box of the object

    When using the GetClosest() method, the user must supply a functor object
    (whose type is a method template argument) which expose the following
    operator ():

    bool operator () (const ObjType & obj, const Point3f & point, ScalarType & mindist, Point3f & result);
    which return true if the distance from point to the object 'obj' is < mindist
    and set mindist to said distance, and result must be set as the closest
    point of the object to point)
    */

    template < class OBJTYPE, class FLT = float >
    class GridStaticPtr : public BasicGrid<FLT>, SpatialIndex<OBJTYPE, FLT>
    {
    public:
        typedef OBJTYPE ObjType;
        typedef ObjType* ObjPtr;
        typedef typename ObjType::ScalarType ScalarType;
        typedef Point3<ScalarType> CoordType;

        typedef GridStaticPtr<OBJTYPE, FLT> GridPtrType;
        typedef BasicGrid<FLT> BT;

        /** Internal class for keeping the first pointer of object.
        Definizione Link dentro la griglia. Classe di supporto per GridStaticObj.
        */
        class Link
        {
        public:
            /// Costruttore di default
            inline Link() {};
            /// Costruttore con inizializzatori
            inline Link(ObjPtr nt, const int ni) {
                assert(ni >= 0);
                t = nt;
                i = ni;
            };


            inline bool operator <  (const Link & l) const { return i < l.i; }
            inline bool operator <= (const Link & l) const { return i <= l.i; }
            inline bool operator >  (const Link & l) const { return i > l.i; }
            inline bool operator >= (const Link & l) const { return i >= l.i; }
            inline bool operator == (const Link & l) const { return i == l.i; }
            inline bool operator != (const Link & l) const { return i != l.i; }

            inline ObjPtr & Elem() {
                return t;
            }

            ObjType &operator *() { return *(t); }

            inline int & Index() {
                return i;
            }

        private:
            /// Puntatore all'elemento T
            ObjPtr t;
            /// Indirizzo del voxel dentro la griglia
            int i;


        };//end class Link

        typedef Link* Cell;
        typedef Cell CellIterator;

        std::vector<Link>   links;   /// Insieme di tutti i links

        std::vector<Cell> grid;   /// Griglia vera e propria

        
        bool Empty() const { return links.empty(); }


        /// BY INTEGER COORDS
        inline Cell* Grid(const int x, const int y, const int z)
        {
            assert(!(x < 0 || x >= BT::siz[0] || y < 0 || y >= BT::siz[1] || z < 0 || z >= BT::siz[2]));
            assert(grid.size() > 0);
            return &*grid.begin() + (x + BT::siz[0] * (y + BT::siz[1] * z));
        }

        void Grid(const int x, const int y, const int z, Cell & first, Cell & last)
        {
            Cell* g = Grid(x, y, z);
            first = *g;
            last = *(g + 1);
        }


        template <class OBJITER>
        inline void Set(const OBJITER & _oBegin, const OBJITER & _oEnd, int _size = 0)
        {
            Box3<FLT> _bbox;
            Box3<FLT> b;
            for (OBJITER i = _oBegin; i != _oEnd; ++i)
            {
                (*i).GetBBox(b);
                _bbox.Add(b);
            }
            if (_size == 0)
                _size = (int)std::distance<OBJITER>(_oBegin, _oEnd);

            ///inflate the bb calculated
            ScalarType infl = _bbox.Diag() / _size;
            _bbox.min -= tessellation::Point3<FLT>(infl, infl, infl);
            _bbox.max += tessellation::Point3<FLT>(infl, infl, infl);

            Set(_oBegin, _oEnd, _bbox, _size);
        }



        // This function automatically compute a reasonable size for the uniform grid such that the number of cells is
        // the same of the nubmer of elements to be inserted in the grid.
        //
        // Note that the bbox must be already 'inflated' so to be sure that no object will fall on the border of the grid.

        template <class OBJITER>
        inline void Set(const OBJITER & _oBegin, const OBJITER & _oEnd, const Box3x &_bbox, int _size = 0)
        {
            if (_size == 0)
                _size = (int)std::distance<OBJITER>(_oBegin, _oEnd);
            Point3<FLT> _dim = _bbox.max - _bbox.min;
            Point3i _siz;
            BestDim(_size, _dim, _siz);

            Set(_oBegin, _oEnd, _bbox, _siz);
        }


        // This is the REAL LOW LEVEL function

        template <class OBJITER>
        inline void Set(const OBJITER & _oBegin, const OBJITER & _oEnd, const Box3x &_bbox, Point3i _siz)
        {
            OBJITER i;

            this->bbox = _bbox;
            this->siz = _siz;

            // find voxel size starting from the provided bbox and grid size. 

            this->dim = this->bbox.max - this->bbox.min;
            this->voxel[0] = this->dim[0] / this->siz[0];
            this->voxel[1] = this->dim[1] / this->siz[1];
            this->voxel[2] = this->dim[2] / this->siz[2];

            // Allocate the grid (add one more for the final sentinel)
            grid.resize(this->siz[0] * this->siz[1] * this->siz[2] + 1);

            // Insert all the objects into the grid
            links.clear();
            for (i = _oBegin; i != _oEnd; ++i)
            {
                Box3x bb;			// Boundig box del tetraedro corrente
                (*i).GetBBox(bb);
                bb.Intersect(this->bbox);
                if (!bb.IsNull())
                {

                    Box3i ib;		// Boundig box in voxels
                    this->BoxToIBox(bb, ib);
                    int x, y, z;
                    for (z = ib.min[2]; z <= ib.max[2]; ++z)
                    {
                        int bz = z*this->siz[1];
                        for (y = ib.min[1]; y <= ib.max[1]; ++y)
                        {
                            int by = (y + bz)*this->siz[0];
                            for (x = ib.min[0]; x <= ib.max[0]; ++x)
                                // Inserire calcolo cella corrente
                                // if( pt->Intersect( ... )
                                links.push_back(Link(&(*i), by + x));
                        }
                    }
                }
            }
            // Push della sentinella
            /*links.push_back( Link((typename ContainerType::iterator)NULL,
            (grid.size()-1)));*/

            links.push_back(Link(NULL, int(grid.size()) - 1));

            // Ordinamento dei links
            sort(links.begin(), links.end());

            // Creazione puntatori ai links
            typename std::vector<Link>::iterator pl;
            unsigned int pg;
            pl = links.begin();
            for (pg = 0; pg < grid.size(); ++pg)
            {
                assert(pl != links.end());
                grid[pg] = &*pl;
                while ((int)pg == pl->Index())
                {
                    ++pl;
                    if (pl == links.end())
                        break;
                }
            }

        }

        template <class OBJPOINTDISTFUNCTOR, class OBJMARKER>
        ObjPtr  GetClosest(OBJPOINTDISTFUNCTOR & _getPointDistance, OBJMARKER & _marker,
            const typename OBJPOINTDISTFUNCTOR::QueryType & _p, const ScalarType & _maxDist, ScalarType & _minDist, CoordType & _closestPt)
        {
            return (tessellation::GridClosest<GridPtrType, OBJPOINTDISTFUNCTOR, OBJMARKER>(*this, _getPointDistance, _marker, _p, _maxDist, _minDist, _closestPt));
        }

    }; //end class GridStaticPtr


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

} // end namespace
