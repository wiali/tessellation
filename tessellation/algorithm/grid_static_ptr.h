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

#include <algorithm/box3.h>
#include <algorithm/grid_closest.h>

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
    class BasicGrid //:public SpatialIndex<SCALARTYPE>
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

                            /*
                            Derives the right values of Dim and voxel starting
                            from the current values of siz and bbox
                            */
        void ComputeDimAndVoxel()
        {
            this->dim = this->bbox.max - this->bbox.min;
            this->voxel[0] = this->dim[0] / this->siz[0];
            this->voxel[1] = this->dim[1] / this->siz[1];
            this->voxel[2] = this->dim[2] / this->siz[2];
        }
        /* Given a 3D point, returns the coordinates of the cell where the point is
        * @param p is a 3D point
        * @return integer coordinates of the cell
        */
        inline Point3i GridP(const Point3<ScalarType> & p) const
        {
            Point3i pi;
            PToIP(p, pi);
            return pi;
        }

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

        /* Given a cell index return the lower corner of the cell
        * @param integer coordinates pi of the cell
        * @return p is a 3D point representing the lower corner of the cell
        */
        template <class OtherScalarType>
        inline void IPiToPf(const Point3i & pi, Point3<OtherScalarType> &p) const
        {
            p[0] = bbox.min[0] + ((OtherScalarType)pi[0])*voxel[0];
            p[1] = bbox.min[1] + ((OtherScalarType)pi[1])*voxel[1];
            p[2] = bbox.min[2] + ((OtherScalarType)pi[2])*voxel[2];
        }

        /* Given a cell index return the corresponding box
        * @param integer coordinates pi of the cell
        * @return b is the corresponding box in <ScalarType> coordinates
        */
        inline void IPiToBox(const Point3i & pi, Box3x & b) const
        {
            CoordType p;
            p[0] = ((ScalarType)pi[0])*voxel[0];
            p[1] = ((ScalarType)pi[1])*voxel[1];
            p[2] = ((ScalarType)pi[2])*voxel[2];
            p += bbox.min;
            b.min = p;
            b.max = (p + voxel);
        }

        /* Given a cell index return the center of the cell itself
        * @param integer coordinates pi of the cell
        * @return b is the corresponding box in <ScalarType> coordinates
        */inline void IPiToBoxCenter(const Point3i & pi, CoordType & c) const
        {
            CoordType p;
            IPiToPf(pi, p);
            c = p + voxel / ScalarType(2.0);
        }

        // Same of IPiToPf but for the case that you just want to transform
        // from a space to the other.
        template <class OtherScalarType>
        void IPfToPf(const Point3<OtherScalarType> & pi, Point3<OtherScalarType> &p) const
        {
            p[0] = ((OtherScalarType)pi[0])*voxel[0] + bbox.min[0];
            p[1] = ((OtherScalarType)pi[1])*voxel[1] + bbox.min[1];
            p[2] = ((OtherScalarType)pi[2])*voxel[2] + bbox.min[2];
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

    template<class scalar_type>
    void BestDim(const Box3<scalar_type> box, const scalar_type voxel_size, Point3i & dim)
    {
        Point3<scalar_type> box_size = box.max - box.min;
        __int64 elem_num = (__int64)(box_size[0] / voxel_size + 0.5) *(__int64)(box_size[1] / voxel_size + 0.5) * (__int64)(box_size[2] / voxel_size + 0.5);
        BestDim(elem_num, box_size, dim);
    }
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
        Method Set.

        Description:
        The Set method initializes the spatial structure.

        Template Parameters:
        OBJITER:  Objects Container's iterator type.

        Method Parameters:
        _oBegin : [IN] begin objects container's iterator
        _oEnd   : [IN] end objects container's iterator

        Return Value:
        None.

        **************************************************************************/
        template <class OBJITER>
        void Set(const OBJITER & _oBegin, const OBJITER & _oEnd) {
            assert(0);      // this is a base interface.
            (void)_oBegin;  // avoid "unreferenced parameter" compiler warning.
            (void)_oEnd;
        }

        /**************************************************************************
        Method Empty.
        Description:
        check if the spatial structure is empty.

        Return Value:
        true if it is empty.
        **************************************************************************/

        bool Empty() {
            assert(0);      // this is a base interface.
            return true;
        }

        /**************************************************************************
        Method GetClosest.

        Description:
        The GetClosest method finds the closest object given a point.
        It also finds the closest point and minimum distance.

        Template Parameters:
        OBJPOINTDISTFUNCTOR : Object-Point distance functor type;
        this type must implement an operator () with signature
        bool operator () (const ObjType & obj, const CoordType & p, ScalarType & d, CoordType & q)
        where:
        obj [IN] is a reference to the current object being tested,
        p [IN] is the query point,
        d [IN/OUT] is in input the reject distance and in output the closest distance,
        q [OUT] is the closest point.
        The operator returns true if the closest distance is less than input reject distance.
        OBJMARKER           : The type of a marker functor.

        Method Parameters:
        _getPointDistance : [IN] Functor for point-distance calculation.
        _marker           : [IN] Functor for marking objects already tested.
        _p                : [IN] The query point.
        _maxDist          : [IN] Maximum reject distance.
        _minDist          : [OUT] Closest distance.
        _closestPt        : [OUT] Closest point.

        Return Value:
        A pointer to the closest object (if any).

        **************************************************************************/
        template <class OBJPOINTDISTFUNCTOR, class OBJMARKER>
        ObjPtr GetClosest(
            OBJPOINTDISTFUNCTOR & _getPointDistance, OBJMARKER & _marker, const CoordType & _p, const ScalarType & _maxDist,
            ScalarType & _minDist, CoordType & _closestPt) {
            assert(0);
            (void)_getPointDistance;
            (void)_marker;
            (void)_p;
            (void)_maxDist;
            (void)_minDist;
            (void)_closestPt;
            return ((ObjPtr)0);
        }

        /**************************************************************************
        Method GetKClosest.

        Description:
        The GetKClosest method finds the K closest object given a point.
        It also finds the closest points and minimum distances.

        Template Parameters:
        OBJPOINTDISTFUNCTOR : Object-Point distance functor type;
        this type must implement an operator () with signature
        bool operator () (const ObjType & obj, const CoordType & p, ScalarType & d, CoordType & q)
        where:
        obj [IN] is a reference to the current object being tested,
        p [IN] is the query point,
        d [IN/OUT] is in input the reject distance and in output the closest distance,
        q [OUT] is the closest point.
        The operator returns true if the closest distance is less than input reject distance.
        OBJMARKER           : The type of a marker functor.
        OBJPTRCONTAINER     : The type of a object pointers container.
        DISTCONTAINER       : The type of a container which, in return, will contain the closest distances.
        POINTCONTAINER      : The type of a container which, in return, will contain the closest points.

        Method Parameters:
        _getPointDistance : [IN] Functor for point-distance calculation.
        _marker           : [IN] Functor for marking objects already tested.
        _k                : [IN] The number of closest objects to search for.
        _p                : [IN] The query point.
        _maxDist          : [IN] Maximum reject distance.
        _objectPtrs       : [OUT] Container which, in return, will contain pointers to the closest objects.
        _distances        : [OUT] Container which, in return, will contain the closest distances.
        _objectPtrs       : [OUT] Container which, in return, will contain the closest points.

        Return Value:
        The number of closest objects found.

        **************************************************************************/
        template <class OBJPOINTDISTFUNCTOR, class OBJMARKER, class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
        unsigned int GetKClosest(
            OBJPOINTDISTFUNCTOR & _getPointDistance, OBJMARKER & _marker, const unsigned int _k, const CoordType & _p, const ScalarType & _maxDist,
            OBJPTRCONTAINER & _objectPtrs, DISTCONTAINER & _distances, POINTCONTAINER & _points) {
            assert(0);
            (void)_getPointDistance;
            (void)_marker;
            (void)_k;
            (void)_p;
            (void)_maxDist;
            (void)_objectPtrs;
            (void)_distances;
            (void)_points;
            return (0);
        }


        /**************************************************************************
        Method GetInSphere.

        Description:
        The GetInSphere method finds all the objects in the specified sphere

        Template Parameters:
        OBJPOINTDISTFUNCTOR : Object-Point distance functor type;
        this type must implement an operator () with signature
        bool operator () (const ObjType & obj, const CoordType & p, ScalarType & d, CoordType & q)
        where:
        obj [IN] is a reference to the current object being tested,
        p [IN] is the query point,
        d [IN/OUT] is in input the reject distance and in output the closest distance,
        q [OUT] is the closest point.
        The operator returns true if the closest distance is less than input reject distance.
        OBJMARKER           : The type of a marker functor.
        OBJPTRCONTAINER     : The type of a object pointers container.
        DISTCONTAINER       : The type of a container which, in return, will contain the closest distances.
        POINTCONTAINER      : The type of a container which, in return, will contain the closest points.

        Method Parameters:
        _getPointDistance : [IN] Functor for point-distance calculation.
        _marker           : [IN] Functor for marking objects already tested.
        _p                : [IN] The query point.
        _r		          : [IN]  The radius of the specified sphere.
        _objectPtrs       : [OUT] Container which, in return, will contain pointers to the in-sphere objects.
        _distances        : [OUT] Container which, in return, will contain the in-sphere distances.
        _objectPtrs       : [OUT] Container which, in return, will contain the in-sphere nearests points for each object.

        Return Value:
        The number of in-sphere objects found.

        **************************************************************************/
        template <class OBJPOINTDISTFUNCTOR, class OBJMARKER, class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
        unsigned int GetInSphere(
            OBJPOINTDISTFUNCTOR & _getPointDistance, OBJMARKER & _marker, const CoordType & _p, const ScalarType & _r, OBJPTRCONTAINER & _objectPtrs, DISTCONTAINER & _distances, POINTCONTAINER & _points) {
            assert(0);
            (void)_getPointDistance;
            (void)_marker;
            (void)_p;
            (void)_r;
            (void)_objectPtrs;
            (void)_distances;
            (void)_points;
            return (0);
        }

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

        /// Date le coordinate di un grid point (corner minx,miy,minz) ritorna le celle che condividono
        /// l'edge cell che parte dal grid point in direzione axis
        inline void Grid(Point3i p, const int axis,
            std::vector<Cell*> & cl)
        {
#ifndef NDEBUG
            if (p[0]<0 || p[0] > BT::siz[0] ||
                p[1]<0 || p[1]> BT::siz[1] ||
                p[2]<0 || p[2]> BT::siz[2])
                assert(0);
            //return NULL;
            else
#endif
                assert(((unsigned int)p[0] + BT::siz[0] * p[1] + BT::siz[1] * p[2]) < grid.size());

            int axis0 = (axis + 1) % 3;
            int axis1 = (axis + 2) % 3;
            int i, j, x, y;
            x = p[axis0];
            y = p[axis1];
            for (i = std::max(x - 1, 0); i <= std::min(x, BT::siz[axis0] - 1); ++i)
                for (j = std::max(y - 1, 0); j <= std::min(y, this->siz[axis1] - 1); ++j) {
                    p[axis0] = i;
                    p[axis1] = j;
                    cl.push_back(Grid(p[0] + BT::siz[0] * (p[1] + BT::siz[1] * p[2])));
                }
        }


        //////////////// 
        // Official access functions
        //////////////// 
        /// BY CELL
        Cell* Grid(const  int i) {
            return &grid[i];
        }

        void Grid(const Cell* g, Cell & first, Cell & last)
        {
            first = *g;
            last = *(g + 1);
        }

        /// BY INTEGER COORDS
        inline Cell* Grid(const int x, const int y, const int z)
        {
            assert(!(x < 0 || x >= BT::siz[0] || y < 0 || y >= BT::siz[1] || z < 0 || z >= BT::siz[2]));
            assert(grid.size() > 0);
            return &*grid.begin() + (x + BT::siz[0] * (y + BT::siz[1] * z));
        }

        inline Cell* Grid(const Point3i &pi)
        {
            return Grid(pi[0], pi[1], pi[2]);
        }

        void Grid(const int x, const int y, const int z, Cell & first, Cell & last)
        {
            Cell* g = Grid(x, y, z);
            first = *g;
            last = *(g + 1);
        }

        void Grid(const Point3<ScalarType> & p, Cell & first, Cell & last)
        {
            Cell* g = Grid(GridP(p));

            first = *g;
            last = *(g + 1);
        }


        /// Set the bounding box of the grid
        ///We need some extra space for numerical precision.
        template <class Box3Type>
        void SetBBox(const Box3Type & b)
        {
            this->bbox.Import(b);
            ScalarType t = this->bbox.Diag() / 100.0;
            if (t == 0) t = ScalarType(1e-20);  // <--- Some doubts on this (Cigno 5/1/04)
            this->bbox.Offset(t);
            this->dim = this->bbox.max - this->bbox.min;
        }



        //		void ShowStats(FILE *fp)
        //		{
        //			// Conto le entry
        //			//int nentry = 0;
        //			//Hist H;
        //			//H.SetRange(0,1000,1000);
        //			//int pg;
        //			//for(pg=0;pg<grid.size()-1;++pg)
        //			//	if( grid[pg]!=grid[pg+1] )
        //			//	{
        //			//		++nentry;
        //			//		H.Add(grid[pg+1]-grid[pg]);
        //			//	}

        //			//	fprintf(fp,"Uniform Grid: %d x %d x %d (%d voxels), %.1f%% full, %d links \nNon empty Cell Occupancy Distribution Avg: %f (%4.0f %4.0f %4.0f) \n",
        //			//	siz[0],siz[1],siz[2],grid.size()-1,
        //			//	double(nentry)*100.0/(grid.size()-1),links.size(),H.Avg(),H.Percentile(.25),H.Percentile(.5),H.Percentile(.75)
        //			//
        //			//);
        //		}

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



        // This function automatically compute a reasonable size for the uniform grid providing the side (radius) of the cell
        //
        // Note that the bbox must be already 'inflated' so to be sure that no object will fall on the border of the grid.

        template <class OBJITER>
        inline void SetWithRadius(const OBJITER & _oBegin, const OBJITER & _oEnd, FLT _cellRadius)
        {
            Box3<FLT> _bbox;
            Box3<FLT> b;
            for (OBJITER i = _oBegin; i != _oEnd; ++i)
            {
                (*i).GetBBox(b);
                _bbox.Add(b);
            }

            _bbox.min -= vcg::Point3<FLT>(_cellRadius, _cellRadius, _cellRadius);
            _bbox.max += vcg::Point3<FLT>(_cellRadius, _cellRadius, _cellRadius);

            Point3i _siz;
            Point3<FLT> _dim = _bbox.max - _bbox.min;
            _dim /= _cellRadius;
            assert(_dim[0] > 0 && _dim[1] > 0 && _dim[2] > 0);
            _siz[0] = (int)ceil(_dim[0]);
            _siz[1] = (int)ceil(_dim[1]);
            _siz[2] = (int)ceil(_dim[2]);

            Set(_oBegin, _oEnd, _bbox, _siz);
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
                while ((int)pg == pl->Index())	// Trovato inizio
                {
                    ++pl;		// Ricerca prossimo blocco
                    if (pl == links.end())
                        break;
                }
            }

        }


        int MemUsed()
        {
            return sizeof(GridStaticPtr) + sizeof(Link)*links.size() +
                sizeof(Cell) * grid.size();
        }

        template <class OBJPOINTDISTFUNCTOR, class OBJMARKER>
        ObjPtr  GetClosest(OBJPOINTDISTFUNCTOR & _getPointDistance, OBJMARKER & _marker,
            const typename OBJPOINTDISTFUNCTOR::QueryType & _p, const ScalarType & _maxDist, ScalarType & _minDist, CoordType & _closestPt)
        {
            return (tessellation::GridClosest<GridPtrType, OBJPOINTDISTFUNCTOR, OBJMARKER>(*this, _getPointDistance, _marker, _p, _maxDist, _minDist, _closestPt));
        }


        template <class OBJPOINTDISTFUNCTOR, class OBJMARKER, class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
        unsigned int GetKClosest(OBJPOINTDISTFUNCTOR & _getPointDistance, OBJMARKER & _marker,
            const unsigned int _k, const CoordType & _p, const ScalarType & _maxDist, OBJPTRCONTAINER & _objectPtrs,
            DISTCONTAINER & _distances, POINTCONTAINER & _points)
        {
            return (vcg::GridGetKClosest<GridPtrType,
                OBJPOINTDISTFUNCTOR, OBJMARKER, OBJPTRCONTAINER, DISTCONTAINER, POINTCONTAINER>(*this, _getPointDistance, _marker, _k, _p, _maxDist, _objectPtrs, _distances, _points));
        }
        
        /* If the grid has a cubic voxel of side <radius> this function
         process all couple of elementes in neighbouring cells.
             GATHERFUNCTOR needs to expose this method:
           bool operator()(OBJTYPE *v1, OBJTYPE *v2);
           which is then called ONCE per unordered pair v1,v2.
         example:

         struct GFunctor {
           double radius2, iradius2;
           GFunctor(double radius) { radius2 = radius*radius; iradius2 = 1/radius2; }

           bool operator()(CVertex *v1, CVertex *v2) {
             Point3d &p = v1->P();
             Point3d &q = v2->P();
             double dist2 = (p-q).SquaredNorm();
             if(dist2 < radius2) {
               double w = exp(dist2*iradius2);
               //do something
             }
           }
        }; */

        template <class GATHERFUNCTOR>
        void Gather(GATHERFUNCTOR gfunctor) {
            static int corner[8 * 3] = { 0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1,
                                       0, 1, 1,  1, 0, 1,  1, 1, 0,  1, 1, 1 };

            static int diagonals[14 * 2] = { 0, 0,
                                           0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7,
                                           2, 3, 1, 3, 1, 2,
                                           1, 4, 2, 5, 3, 6 };

            Cell ostart, oend, dstart, dend;
            for (int z = 0; z < this->siz[2]; z++) {
                for (int y = 0; y < this->siz[1]; y++) {
                    for (int x = 0; x < this->siz[0]; x++) {

                        Grid(x, y, z, ostart, oend);

                        for (Cell c = ostart; c != oend; c++)
                            for (Cell s = c + 1; s != oend; s++)
                                gfunctor(c->Elem(), s->Elem());

                        for (int d = 2; d < 28; d += 2) { //skipping self
                            int *cs = corner + 3 * diagonals[d];
                            int *ce = corner + 3 * diagonals[d + 1];
                            if ((x + cs[0] < this->siz[0]) && (y + cs[1] < this->siz[1]) && (z + cs[2] < this->siz[2]) &&
                                (x + ce[0] < this->siz[0]) && (y + ce[1] < this->siz[1]) && (z + ce[2] < this->siz[2])) {

                                Grid(x + cs[0], y + cs[1], z + cs[2], ostart, oend);
                                Grid(x + ce[0], y + ce[1], z + ce[2], dstart, dend);

                                for (Cell c = ostart; c != oend; c++)
                                    for (Cell s = dstart; s != dend; s++)
                                        gfunctor(c->Elem(), s->Elem());
                            }
                        }
                    }
                }
            }
        }


    }; //end class GridStaticPtr

} // end namespace
