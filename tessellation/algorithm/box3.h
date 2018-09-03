#pragma once

#include <algorithm/space_base.h>

namespace tessellation
{

    template <class BoxScalarType>
    class Box3
    {
    public:

        /// The scalar type
        typedef BoxScalarType ScalarType;

        /// min coordinate point
        Point3<BoxScalarType> min;
        /// max coordinate point
        Point3<BoxScalarType> max;
        /// The bounding box constructor
        inline  Box3() { min.X() = 1; max.X() = -1; min.Y() = 1; max.Y() = -1; min.Z() = 1; max.Z() = -1; }
        /// Copy constructor
        inline  Box3(const Box3 & b) { min = b.min; max = b.max; }
        /// Min Max constructor
        inline  Box3(const Point3<BoxScalarType> & mi, const Point3<BoxScalarType> & ma) { min = mi; max = ma; }
        /// Point Radius Constructor
        inline Box3(const Point3<BoxScalarType> & center, const BoxScalarType & radius) {
            min = center - Point3<BoxScalarType>(radius, radius, radius);
            max = center + Point3<BoxScalarType>(radius, radius, radius);
        }
        /// The bounding box distructor
        inline ~Box3() { }
        /// Operator to compare two bounding box
        inline bool operator == (Box3<BoxScalarType> const & p) const
        {
            return min == p.min && max == p.max;
        }
        /// Operator to dispare two bounding box
        inline bool operator != (Box3<BoxScalarType> const & p) const
        {
            return min != p.min || max != p.max;
        }
        /** Offset of a vector (s,s,s)
        */
        void Offset(const BoxScalarType s)
        {
            Offset(Point3<BoxScalarType>(s, s, s));
        }
        /** Offset the two corner of the box of a vector delta.
         *  adding delta to max and -delta to min.
            @param delta offset vector
        */
        void Offset(const Point3<BoxScalarType> & delta)
        {
            min -= delta;
            max += delta;
        }
        /// Initializing the bounding box
        void Set(const Point3<BoxScalarType> & p)
        {
            min = max = p;
        }

        /// Set the bounding box to a null value
        void SetNull()
        {
            min.X() = 1; max.X() = -1;
            min.Y() = 1; max.Y() = -1;
            min.Z() = 1; max.Z() = -1;
        }
        /** Modify the current bbox to contain also the passed box.
         *  Adding a null bounding box does nothing
        */
        void Add(Box3<BoxScalarType> const & b)
        {
            if (b.IsNull()) return; // Adding a null bbox should do nothing
            if (IsNull()) *this = b;
            else
            {
                if (min.X() > b.min.X()) min.X() = b.min.X();
                if (min.Y() > b.min.Y()) min.Y() = b.min.Y();
                if (min.Z() > b.min.Z()) min.Z() = b.min.Z();

                if (max.X() < b.max.X()) max.X() = b.max.X();
                if (max.Y() < b.max.Y()) max.Y() = b.max.Y();
                if (max.Z() < b.max.Z()) max.Z() = b.max.Z();
            }
        }
        /** Modify the current bbox to contain also the passed point
        */
        void Add(const Point3<BoxScalarType> & p)
        {
            if (IsNull()) Set(p);
            else
            {
                if (min.X() > p.X()) min.X() = p.X();
                if (min.Y() > p.Y()) min.Y() = p.Y();
                if (min.Z() > p.Z()) min.Z() = p.Z();

                if (max.X() < p.X()) max.X() = p.X();
                if (max.Y() < p.Y()) max.Y() = p.Y();
                if (max.Z() < p.Z()) max.Z() = p.Z();
            }
        }

        /** Modify the current bbox to contain also the passed sphere
        */
        void Add(const Point3<BoxScalarType> & p, const BoxScalarType radius)
        {
            if (IsNull()) Set(p);
            else
            {
                min.X() = std::min(min.X(), p.X() - radius);
                min.Y() = std::min(min.Y(), p.Y() - radius);
                min.Z() = std::min(min.Z(), p.Z() - radius);

                max.X() = std::max(max.X(), p.X() + radius);
                max.Y() = std::max(max.Y(), p.Y() + radius);
                max.Z() = std::max(max.Z(), p.Z() + radius);
            }
        }

        void Intersect(const Box3<BoxScalarType> & b)
        {
            if (min.X() < b.min.X()) min.X() = b.min.X();
            if (min.Y() < b.min.Y()) min.Y() = b.min.Y();
            if (min.Z() < b.min.Z()) min.Z() = b.min.Z();

            if (max.X() > b.max.X()) max.X() = b.max.X();
            if (max.Y() > b.max.Y()) max.Y() = b.max.Y();
            if (max.Z() > b.max.Z()) max.Z() = b.max.Z();

            if (min.X() > max.X() || min.Y() > max.Y() || min.Z() > max.Z()) SetNull();
        }
        /** Trasla il bounding box di un valore definito dal parametro.
            @param p Il bounding box trasla sulla x e sulla y in base alle coordinate del parametro
        */
        void Translate(const Point3<BoxScalarType> & p)
        {
            min += p;
            max += p;
        }
        /** true if the point belong to the closed box
        */
        bool IsIn(Point3<BoxScalarType> const & p) const
        {
            return (
                min.X() <= p.X() && p.X() <= max.X() &&
                min.Y() <= p.Y() && p.Y() <= max.Y() &&
                min.Z() <= p.Z() && p.Z() <= max.Z()
                );
        }
        /** true if the point belong to the open box (open on the max side)
         * e.g. if p in [min,max)
        */
        bool IsInEx(Point3<BoxScalarType> const & p) const
        {
            return (
                min.X() <= p.X() && p.X() < max.X() &&
                min.Y() <= p.Y() && p.Y() < max.Y() &&
                min.Z() <= p.Z() && p.Z() < max.Z()
                );
        }
        /** Verifica se due bounding box collidono cioe' se hanno una intersezione non vuota. Per esempio
            due bounding box adiacenti non collidono.
            @param b A bounding box
            @return True se collidoo, false altrimenti
        */
        /* old version
        bool Collide(Box3<BoxScalarType> const &b)
        {
            Box3<BoxScalarType> bb=*this;
            bb.Intersect(b);
            return bb.IsValid();
        }
        */
        bool Collide(Box3<BoxScalarType> const &b) const
        {
            return b.min.X() < max.X() && b.max.X() > min.X() &&
                b.min.Y() < max.Y() && b.max.Y() > min.Y() &&
                b.min.Z() < max.Z() && b.max.Z() > min.Z();
        }
        /**
          return true if the box is null (e.g. invalid or not initialized);
        */
        bool IsNull() const { return min.X() > max.X() || min.Y() > max.Y() || min.Z() > max.Z(); }
        /** return true if the box is empty (e.g. if min == max)
        */
        bool IsEmpty() const { return min == max; }
        /// Return the lenght of the diagonal of the box .
        BoxScalarType Diag() const
        {
            return Distance(min, max);
        }
        /// Calcola il quadrato della diagonale del bounding box.
        BoxScalarType SquaredDiag() const
        {
            return SquaredDistance(min, max);
        }
        /// Return the center of the box.
        Point3<BoxScalarType> Center() const
        {
            return (min + max) / 2;
        }
        /// Compute bounding box size.
        Point3<BoxScalarType> Dim() const
        {
            return (max - min);
        }
        /// Returns global coords of a local point expressed in [0..1]^3
        Point3<BoxScalarType> LocalToGlobal(Point3<BoxScalarType> const & p) const {
            return Point3<BoxScalarType>(
                min[0] + p[0] * (max[0] - min[0]),
                min[1] + p[1] * (max[1] - min[1]),
                min[2] + p[2] * (max[2] - min[2]));
        }
        /// Returns local coords expressed in [0..1]^3 of a point in 3D
        Point3<BoxScalarType> GlobalToLocal(Point3<BoxScalarType> const & p) const {
            return Point3<BoxScalarType>(
                (p[0] - min[0]) / (max[0] - min[0]),
                (p[1] - min[1]) / (max[1] - min[1]),
                (p[2] - min[2]) / (max[2] - min[2])
                );
        }
        /// Return the volume of the box.
        BoxScalarType Volume() const
        {
            return (max.X() - min.X())*(max.Y() - min.Y())*(max.Z() - min.Z());
        }
        /// Calcola la dimensione del bounding box sulla x.
        inline BoxScalarType DimX() const { return max.X() - min.X(); }
        /// Calcola la dimensione del bounding box sulla y.
        inline BoxScalarType DimY() const { return max.Y() - min.Y(); }
        /// Calcola la dimensione del bounding box sulla z.
        inline BoxScalarType DimZ() const { return max.Z() - min.Z(); }
        /// Calcola il lato di lunghezza maggiore
        inline unsigned char MaxDim() const {
            int i;
            Point3<BoxScalarType> diag = max - min;
            if (diag[0] > diag[1]) i = 0; else i = 1;
            return (diag[i] > diag[2]) ? i : 2;
        }
        /// Calcola il lato di lunghezza minore
        inline unsigned char MinDim() const {
            int i;
            Point3<BoxScalarType> diag = max - min;
            if (diag[0] < diag[1]) i = 0; else i = 1;
            return (diag[i] < diag[2]) ? i : 2;
        }

        template <class Q>
        inline void Import(const Box3<Q> & b)
        {
            min.Import(b.min);
            max.Import(b.max);
        }

        template <class Q>
        static inline Box3 Construct(const Box3<Q> & b)
        {
            return Box3(Point3<BoxScalarType>::Construct(b.min), Point3<BoxScalarType>::Construct(b.max));
        }

        /// gives the ith box vertex in order: (x,y,z),(X,y,z),(x,Y,z),(X,Y,z),(x,y,Z),(X,y,Z),(x,Y,Z),(X,Y,Z)
        Point3<BoxScalarType> P(const int & i) const {
            return Point3<BoxScalarType>(
                min[0] + (i % 2) * DimX(),
                min[1] + ((i / 2) % 2) * DimY(),
                min[2] + (i > 3)* DimZ());
        }
    }; // end class definition



    typedef Box3<short>  Box3s;
    typedef Box3<int>	 Box3i;
    typedef Box3<float>  Box3f;
    typedef Box3<double> Box3d;


    template <class LineScalarType, bool NORM = false>
    class Line3
    {
    public:

        /// The scalar type
        typedef LineScalarType ScalarType;

        /// The point type
        typedef Point3<LineScalarType> PointType;

        /// The line type
        typedef Line3<LineScalarType, NORM> LineType;

    private:

        /// Origin
        PointType _ori;

        /// Direction (not necessarily normalized, unless so specified by NORM)
        PointType _dir;

    public:

        //@{
        /** @name Members to access the origin or direction
        Direction() cannot be assigned directly.
        Use SetDirection() or Set() instead.
        **/
        /// 
        inline const PointType &Origin() const { return _ori; }
        inline PointType &Origin() { return _ori; }
        inline const PointType &Direction() const { return _dir; }
        /// sets the origin
        inline void SetOrigin(const PointType & ori)
        {
            _ori = ori;
        }
        /// sets the direction
        inline void SetDirection(const PointType & dir)
        {
            _dir = dir; if (NORM) _dir.Normalize();
        }
        /// sets origin and direction.
        inline void Set(const PointType & ori, const PointType & dir)
        {
            SetOrigin(ori); SetDirection(dir);
        }
        //@}

        //@{
        /** @name Constructors
        **/
        /// The empty constructor
        Line3() {};
        /// The (origin, direction) constructor
        Line3(const PointType &ori, const PointType &dir) { SetOrigin(ori); SetDirection(dir); };
        //@}

        /// Operator to compare two lines
        inline bool operator == (LineType const & p) const
        {
            return _ori == p._ori && _dir == p._dir;
        }
        /// Operator to dispare two lines
        inline bool operator != (LineType const & p) const
        {
            return _ori != p._ori || _dir != p._dir;
        }
        /// Projects a point on the line
        inline ScalarType Projection(const  PointType &p) const
        {
            if (NORM) return ScalarType((p - _ori).dot(_dir));
            else      return ScalarType((p - _ori).dot(_dir) / _dir.SquaredNorm());
        }
        /// returns wheter this type is normalized or not
        static bool IsNormalized() { return NORM; };
        /// calculates the point of parameter t on the line.
        inline PointType P(const ScalarType t) const
        {
            return _ori + _dir * t;
        }
        /// normalizes direction field (returns a Normalized Line)
        inline Line3<ScalarType, true> &Normalize()
        {
            if (!NORM) _dir.Normalize(); return *((Line3<ScalarType, true>*)this);
        }
        /// normalizes direction field (returns a Normalized Line) - static version
        static Line3<ScalarType, true> &Normalize(LineType &p)
        {
            p.Normalize(); return *((Line3<ScalarType, true>*)(&p));
        }
        /// importer for different line types (with any scalar type or normalization beaviour)
        template <class Q, bool K>
        inline void Import(const Line3<Q, K> & b)
        {
            _ori.Import(b.Origin());	_dir.Import(b.Direction());
            if ((NORM) && (!K)) _dir.Normalize();
            //printf("(=)%c->%c ",(!NORM)?'N':'n', NORM?'N':'n');
        }
        /// constructs a new line importing it from an existing one
        template <class Q, bool K>
        static LineType Construct(const Line3<Q, K> & b)
        {
            LineType res; res.Import(b);  return res;
        }
        PointType ClosestPoint(const PointType & p) const {
            return P(Projection(p));
        }
        /// flips the line
        inline void Flip() {
            _dir = -_dir;
        };

        //@{
        /** @name Linearity for 3d lines
        (operators +, -, *, /) so a line can be set as a linear combination
        of several lines. Note that the result of any operation returns
        a non-normalized line; however, the command r0 = r1*a + r2*b is licit
        even if r0,r1,r2 are normalized lines, as the normalization will
        take place within the final assignement operation.
        **/
        inline Line3<ScalarType, false> operator + (LineType const & p) const
        {
            return Line3<ScalarType, false>(_ori + p.Origin(), _dir + p.Direction());
        }
        inline Line3<ScalarType, false> operator - (LineType const & p) const
        {
            return Line3<ScalarType, false>(_ori - p.Origin(), _dir - p.Direction());
        }
        inline Line3<ScalarType, false> operator * (const ScalarType s) const
        {
            return Line3<ScalarType, false>(_ori*s, _dir*s);
        }
        inline Line3<ScalarType, false> operator / (const ScalarType s) const
        {
            ScalarType s0 = ((ScalarType)1.0) / s; return LineType(_ori*s0, _dir*s0);
        }
        //@}


        //@{
        /** @name Automatic normalized to non-normalized
        "Line3dN r0 = r1" is equivalent to
        "Line3dN r0 = r1.Normalize()" if r1 is a Line3d
        **/
        /// copy constructor that takes opposite beaviour
        Line3(const Line3<ScalarType, !NORM > &r)
        {
            Import(r);
        };
        /// assignment
        inline LineType & operator = (Line3<ScalarType, !NORM> const &r)
        {
            Import(r); return *this;
        };
        //@}

    }; // end class definition

    typedef Line3<short>  Line3s;
    typedef Line3<int>	  Line3i;
    typedef Line3<float>  Line3f;
    typedef Line3<double> Line3d;

    typedef Line3<short, true> Line3sN;
    typedef Line3<int, true> Line3iN;
    typedef Line3<float, true> Line3fN;
    typedef Line3<double, true> Line3dN;

    /// returns closest point
    template <class ScalarType, bool NORM>
    Point3<ScalarType> ClosestPoint(Line3<ScalarType, NORM> l, const Point3<ScalarType> & p)
    {
        return l.P(l.Projection(p));
    }

    template <class ScalarType, bool NORM>
    ScalarType Distance(const Line3<ScalarType, NORM> &l,
        const Point3<ScalarType> &p) {
        Point3<ScalarType> o = l.ClosestPoint(p);
        return (o - p).Norm();
    }

    /*@}*/

    /** \addtogroup space */
    /*@{*/
    /**
    Templated class for 3D segment.
    This is the class for a segment in 3D space. A Segment is stored just as its two extrema (Point3).
    @param SegmentScalarType (template parameter) Specifies the type of scalar used to represent coords.
    */
    template <class SegmentScalarType >
    class Segment3
    {
    public:

        /// The scalar type
        typedef SegmentScalarType ScalarType;

        /// The point type
        typedef Point3<SegmentScalarType> PointType;

        /// The point type
        typedef Segment3<SegmentScalarType> SegmentType;

    private:

        /// _extrema
        PointType _p0, _p1;

    public:

        /// Members to access either extrema
        inline const PointType &P0() const { return _p0; }
        inline const PointType &P1() const { return _p1; }
        inline PointType &P0() { return _p0; }
        inline PointType &P1() { return _p1; }

        inline PointType &operator[] (const int i) { return i == 0 ? _p0 : _p1; }
        inline const PointType &operator[] (const int i) const { return i == 0 ? _p0 : _p1; }

        /// The empty constructor
        Segment3() {};
        /// The (a,b) constructor
        Segment3(const PointType &a, const PointType &b) { _p0 = a; _p1 = b; };
        /// Operator to compare segments
        inline bool operator == (SegmentType const & p) const
        {
            return _p0 == p._p0 && _p1 == p._p1;
        }
        /// Operator to dispare segments
        inline bool operator != (SegmentType const & p) const
        {
            return _p0 != p._p0 || _p1 != p._p1;
        }
        /// initializes the segment with its extrema
        void Set(const PointType &a, const PointType &b)
        {
            _p0 = a; _p1 = b;
        }
        /// calculates the point of parameter t on the segment.
        /// if t is in [0..1] returned point is inside the segment
        inline PointType Lerp(const ScalarType t) const
        {
            return _p0 + (_p1 - _p0) * t;
        }
        /// return the middle point
        inline PointType MidPoint() const
        {
            return (_p0 + _p1) / ScalarType(2.0);
        }
        inline PointType Direction() const
        {
            return (_p1 - _p0);
        }
        inline PointType NormalizedDirection() const
        {
            return (_p1 - _p0).Normalize();
        }
        /// return the bounding box
        inline Box3<ScalarType> BBox() const
        {
            Box3<ScalarType> t;
            if (_p0[0]<_p1[0]) { t.min[0] = _p0[0]; t.max[0] = _p1[0]; }
            else { t.min[0] = _p1[0]; t.max[0] = _p0[0]; }
            if (_p0[1]<_p1[1]) { t.min[1] = _p0[1]; t.max[1] = _p1[1]; }
            else { t.min[1] = _p1[1]; t.max[1] = _p0[1]; }
            if (_p0[2]<_p1[2]) { t.min[2] = _p0[2]; t.max[2] = _p1[2]; }
            else { t.min[2] = _p1[2]; t.max[2] = _p0[2]; }
            return t;
        }
        /// returns segment length
        ScalarType Length() const
        {
            return (_p0 - _p1).Norm();
        }
        /// return segment squared length
        ScalarType SquaredLength() const
        {
            return (_p0 - _p1).SquaredNorm();
        }
        /// flips: a-b becomes b-a
        void Flip()
        {
            PointType t = _p0; _p0 = _p1; _p1 = t;
        }
        /// importer for different line types
        template <class Q>
        inline void Import(const Segment3<Q> & b)
        {
            _p0.Import(b.P0());	_p1.Import(b.P1());
        }
        /// copy constructor (builds a new segment importing an existing one)
        template <class Q>
        static SegmentType Construct(const Segment3<Q> & b)
        {
            return SegmentType(PointType::Construct(b.P0()), PointType::Construct(b.P1()));
        }

        //@{
        /** @name Linearity for 3d segments (operators +, -, *, /) **/
        inline SegmentType operator + (SegmentType const & p) const
        {
            return SegmentType(_p0 + p.P0(), _p1 + p.P1());
        }
        inline SegmentType operator - (SegmentType const & p) const
        {
            return SegmentType(_p0 - p.P0(), _p1 - p.P1());
        }
        inline SegmentType operator * (const ScalarType s) const
        {
            return SegmentType(_p0*s, _p1*s);
        }
        inline SegmentType operator / (const ScalarType s) const
        {
            ScalarType s0 = ((ScalarType)1.0) / s; return SegmentType(_p0*s0, _p1*s0);
        }
        //@}

    }; // end class definition

    typedef Segment3<short>  Segment3s;
    typedef Segment3<int>	   Segment3i;
    typedef Segment3<float>  Segment3f;
    typedef Segment3<double> Segment3d;


    /*
    * Computes the minimum distance between a segment and a point
    * @param[in] segment	The input segment
    * @param[in] p				The input point
    * @param[in] clos			The closest point
    * @param[in] dist The distance
    */
    template <class ScalarType>
    void SegmentPointDistance(Segment3<ScalarType> s,
        const Point3<ScalarType> & p,
        Point3< ScalarType > &clos,
        ScalarType &dist)
    {
        SegmentPointSquaredDistance(s, p, clos, dist);
        dist = sqrt(dist);
    }


    /*
    * Computes the minimum distance between a segment and a point
    * @param[in] segment	The input segment
    * @param[in] p				The input point
    * @param[in] clos			The closest point
    * @param[in] sqr_dist The squared distance
    */
    template <class ScalarType>
    void SegmentPointSquaredDistance(const Segment3<ScalarType> &s,
        const Point3<ScalarType> & p,
        Point3< ScalarType > &closest,
        ScalarType &sqr_dist)
    {
        Point3<ScalarType> e = s.P1() - s.P0();
        ScalarType eSquaredNorm = e.SquaredNorm();
        if (eSquaredNorm < std::numeric_limits<ScalarType>::min())
        {
            closest = s.MidPoint();
            sqr_dist = SquaredDistance(closest, p);
        }
        else
        {
            ScalarType  t = ((p - s.P0())*e) / eSquaredNorm;
            if (t < 0)      t = 0;
            else if (t > 1) t = 1;
            closest = s.P0() + e*t;
            sqr_dist = SquaredDistance(p, closest);
            assert(!math::IsNAN(sqr_dist));
        }
    }

    template <class SCALAR_TYPE> class Triangle2
    {
    public:
        typedef SCALAR_TYPE ScalarType;
        typedef Point2< ScalarType > CoordType;
        typedef Triangle2<ScalarType> TriangleType;

    protected:
        /// Vector of vertex pointer incident in the face
        Point2<ScalarType> _v[3];
    public:

        Triangle2()
        {}

        Triangle2(const CoordType &p0, const CoordType &p1, const CoordType &p2)
        {
            P(0) = p0;
            P(1) = p1;
            P(2) = p2;
        }

        /// Shortcut per accedere ai punti delle facce
        inline CoordType & P(const int j) { return _v[j]; }
        inline CoordType & P0(const int j) { return _v[j]; }
        inline CoordType & P1(const int j) { return _v[(j + 1) % 3]; }
        inline CoordType & P2(const int j) { return _v[(j + 2) % 3]; }
        inline const CoordType &  P(const int j) const { return _v[j]; }
        inline const CoordType &  P0(const int j) const { return _v[j]; }
        inline const CoordType &  P1(const int j) const { return _v[(j + 1) % 3]; }
        inline const CoordType &  P2(const int j) const { return _v[(j + 2) % 3]; }
        inline const CoordType & cP0(const int j) const { return _v[j]; }
        inline const CoordType & cP1(const int j) const { return _v[(j + 1) % 3]; }
        inline const CoordType & cP2(const int j) const { return _v[(j + 2) % 3]; }

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
        bool InterpolationParameters(const CoordType & bq, ScalarType &L1,
            ScalarType &L2, ScalarType &L3) const
        {
            const ScalarType EPSILON = ScalarType(0.0001f);

            ScalarType x1 = P(0).X();
            ScalarType x2 = P(1).X();
            ScalarType x3 = P(2).X();

            ScalarType y1 = P(0).Y();
            ScalarType y2 = P(1).Y();
            ScalarType y3 = P(2).Y();

            ScalarType x = bq.X();
            ScalarType y = bq.Y();

            L1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
            L2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / ((y3 - y1)*(x2 - x3) + (x1 - x3)*(y2 - y3));
            L3 = 1 - L1 - L2;
            if (math::IsNAN(L1) || math::IsNAN(L2) || math::IsNAN(L3)) L1 = L2 = L3 = (ScalarType)(1.0 / 3.0);
            bool inside = true;
            inside &= (L1 >= 0 - EPSILON) && (L1 <= 1 + EPSILON);
            inside &= (L2 >= 0 - EPSILON) && (L2 <= 1 + EPSILON);
            inside &= (L3 >= 0 - EPSILON) && (L3 <= 1 + EPSILON);
            return inside;
        }

        ///return the distance to the point q and neighors point p
        void PointDistance(const CoordType & q,
            ScalarType & dist,
            CoordType & p) const
        {
            dist = FLT_MAX;
            ///find distance to each segment and take minimum
            for (int i = 0; i<3; i++)
            {
                vcg::Segment2<ScalarType> s = vcg::Segment2<ScalarType>(P(i), P((i + 1) % 3));
                CoordType clos = ClosestPoint<ScalarType>(s, q);
                ScalarType dis_test = (clos - q).Norm();
                if (dis_test<dist)
                {
                    dist = dis_test;
                    p = clos;
                }
            }
        }

        ///retutn true if the face is contuerclockwise oriented
        bool IsCCW()
        {
            ScalarType Area = (P(1) - P(0)) ^ (P(2) - P(0));
            return (Area>0);
        }

    }; //end Class

       /** \addtogroup space */
       /*@{*/
       /**
       Templated class for 2D planes in 3D spaces.
       This is the class for infinite planes in 3D space. A Plane is stored just as a Point3 and a scalar:
       * a direction (not necessarily normalized),
       * an offset from the origin

       Just to be clear, given a point P on a plane it always holds:

       plane.Direction().dot(P) == plane.Offset()


       @param T (template parameter) Specifies the type of scalar used to represent coords.
       @param NORM: if on, the direction is always Normalized
       */
    template <class T, bool NORM = true> class Plane3 {
    public:
        typedef T ScalarType;
        typedef Point3<T> PointType;

    private:
        /// Distance
        ScalarType _offset;
        ///Direction (not necessarily normalized unless NORM is true)
        PointType _dir;

    public:
        //@{
        /** @name Constructors
        **/
        /// The empty constructor
        Plane3() {}
        /// The (distance, direction) constructor
        Plane3(const ScalarType &dist, const PointType &dir) { Set(dist, dir); }

        template <class Q>
        inline void Import(const Plane3<Q, false> & b)
        {
            _offset = ScalarType(b.Offset());
            _dir = Point3<T>::Construct(b.Direction());
        }

        //@{
        /** @name Members to access the distance or direction
        Direction() cannot be assigned directly.
        Use SetDirection() or Set() instead. This is mandatory to make possible the automatic autonormalization template mechanism.
        Note that if you have to set both direction and offset it can be more efficient to set them toghether
        **/
        const ScalarType &Offset() const { return _offset; }
        ScalarType &Offset() { return _offset; }
        /// sets the origin
        void SetOffset(const ScalarType &o) { _offset = o; }

        const PointType &Direction() const { return _dir; }
        /// sets the direction
        void SetDirection(const PointType & dir) {
            _dir = dir;
            if (NORM) _dir.Normalize();
        }
        /// sets origin and direction.
        void Set(const ScalarType & off, const PointType & dir) {
            if (NORM) {
                const ScalarType normFactor = dir.Norm();
                this->_dir = dir / normFactor;
                this->_offset = off / normFactor;
            }
            else {
                this->_offset = off;
                this->_dir = dir;
            }
        }
        void Set(const PointType & dir, const ScalarType & off) { Set(off, dir); }

        /// Operator to compare two lines
        bool operator==(Plane3 const &p) const {
            return _offset == p._offset && _dir == p._dir;
        }
        /// Operator to dispare two lines
        bool operator!=(Plane3 const &p) const {
            return _offset != p._offset || _dir != p._dir;
        }

        ///Project a point on the plane
        PointType Projection(const PointType &p) const {
            ScalarType k = p.dot(_dir) - _offset;
            return p - _dir * k;
        }

        ///Mirror the point wrt the plane
        PointType Mirror(const PointType &p) const {
            PointType mirr = Projection(p);
            mirr += mirr - p;
            return mirr;
        }

        /// Function to normalize direction
        void Normalize() {
            _dir.Normalize();
        }

        /// Calculates the plane passing through three points (Rename this method)
        void Init(const PointType &p0, const PointType &p1, const PointType &p2) {
            _dir = (p2 - p0) ^ (p1 - p0);
            if (NORM) Normalize();
            _offset = p0.dot(_dir);
        }

        /// Calculates the plane passing through a point and the normal (Rename this method
        inline void Init(const PointType &p0, const PointType &norm) {
            _dir = norm;
            if (NORM) Normalize();
            _offset = p0.dot(_dir);
        }
    };	// end class Plane3

    typedef Plane3<float>  Plane3f;
    typedef Plane3<double> Plane3d;

    ///Distance plane - point and vv. (Move these function to somewhere else)
    template<class T> T SignedDistancePlanePoint(const Plane3<T, true> & plane, const Point3<T> & point)
    {
        return plane.Direction().dot(point) - plane.Offset();
    }
    
    template<class T> T SignedDistancePointPlane(const Point3<T> & point, const Plane3<T, true> & plane)
    {
        return SignedDistancePlanePoint(plane, point);
    }
} // end namespace
