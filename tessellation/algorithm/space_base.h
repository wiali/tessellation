#pragma once

#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace tessellation
{
    namespace math
    {
        inline float Sqrt(const float v) { return sqrtf(v); }

        template<class T> int IsNAN(T t) { return _isnan(t) || (!_finite(t)); }

        template<class T> inline const T & Min(const T &a, const T &b, const T &c) {
            if (a < b) {
                if (a < c) return a;
                else return c;
            }
            else {
                if (b < c) return b;
                else return c;
            }
        }

        template<class T> inline const T & Max(const T &a, const T &b, const T &c) {
            if (a > b) {
                if (a > c) return a;
                else return c;
            }
            else {
                if (b > c) return b;
                else return c;
            }
        }
    }

    template < typename T = int>
    class DefaultDeriver : public T {};

    template <
        class Base,
        template <typename> class A>
    class Arity1 : public A<Base > {
    };

    template <
        class Base,
        template <typename> class A, template <typename> class B>
    class Arity2 : public B<Arity1<Base, A> > {};

    template <
        class Base,
        template <typename> class A, template <typename> class B,
        template <typename> class C >
    class Arity3 : public C<Arity2<Base, A, B> > {};

    template <
        class Base,
        template <typename> class A, template <typename> class B,
        template <typename> class C, template <typename> class D>
    class Arity4 : public D<Arity3<Base, A, B, C> > {};

    template <
        class Base,
        template <typename> class A, template <typename> class B,
        template <typename> class C, template <typename> class D,
        template <typename> class E >
    class Arity5 : public E<Arity4<Base, A, B, C, D> > {};

    template <
        class Base,
        template <typename> class A, template <typename> class B,
        template <typename> class C, template <typename> class D,
        template <typename> class E, template <typename> class F >
    class Arity6 : public F<Arity5<Base, A, B, C, D, E> > {};

    template <
        class Base,
        template <typename> class A, template <typename> class B,
        template <typename> class C, template <typename> class D,
        template <typename> class E, template <typename> class F,
        template <typename> class G>
    class Arity7 : public G<Arity6<Base, A, B, C, D, E, F> > {};

    template <
    class Base,
          class TA,
          template <typename, typename> class A >
    class MArity1 : public A<  Base, TA  >
    {
    };

    template <
    class Base,
          class TA,
          template <typename, typename> class A,
          class TB,
          template <typename, typename> class B >
    class MArity2 : public B< MArity1<Base, TA, A>, TB > {};
    template <
        class Base,
        class TA,
        template <typename, typename> class A,
        class TB,
        template <typename, typename> class B,
        class TC,
        template <typename, typename> class C >
    class MArity3 : public C<MArity2<Base, TA, A, TB, B>, TC > {};

    template <
        class Base,
        class TA,
        template <typename, typename> class A,
        class TB,
        template <typename, typename> class B,
        class TC,
        template <typename, typename> class C,
        class TD,
        template <typename, typename> class D>
    class MArity4 : public D<MArity3<Base, TA, A, TB, B, TC, C>, TD > {};

    class DumClass {};

    template <class P2ScalarType> class Point2
    {
    protected:
        P2ScalarType _v[2];
    public:
        typedef P2ScalarType ScalarType;
        enum { Dimension = 2 };

        inline Point2() { }
        inline Point2(const ScalarType nx, const ScalarType ny)
        {
            _v[0] = nx;
            _v[1] = ny;
        }

        inline const ScalarType &X() const { return _v[0]; }
        inline const ScalarType &Y() const { return _v[1]; }
        inline ScalarType &X() { return _v[0]; }
        inline ScalarType &Y() { return _v[1]; }

        inline Point2 operator + (Point2 const &p) const
        {
            return Point2<ScalarType>(_v[0] + p._v[0], _v[1] + p._v[1]);
        }
        inline Point2 operator - (Point2 const &p) const
        {
            return Point2<ScalarType>(_v[0] - p._v[0], _v[1] - p._v[1]);
        }
        inline Point2 operator * (const ScalarType s) const
        {
            return Point2<ScalarType>(_v[0] * s, _v[1] * s);
        }
        inline Point2 operator / (const ScalarType s) const
        {
            return Point2<ScalarType>(_v[0] / s, _v[1] / s);
        }

        inline const ScalarType & operator [] (const int i) const
        {
            assert(i >= 0 && i < 2);
            return _v[i];
        }
        inline ScalarType & operator [] (const int i)
        {
            assert(i >= 0 && i < 2);
            return _v[i];
        }

        inline Point2 & Normalize(void)
        {
            ScalarType n = math::Sqrt(_v[0] * _v[0] + _v[1] * _v[1]);
            if (n > 0.0) { _v[0] /= n;	_v[1] /= n; }
            return *this;
        }

    };

    typedef Point2<float>  Point2f;

    template <class T> class Box3;

    template <class P3ScalarType> class Point3
    {
    protected:
        P3ScalarType _v[3];

    public:
        typedef P3ScalarType ScalarType;
        enum { Dimension = 3 };

        inline Point3() { }
        inline Point3(const P3ScalarType nx, const P3ScalarType ny, const P3ScalarType nz)
        {
            _v[0] = nx;
            _v[1] = ny;
            _v[2] = nz;
        }

        template <class Q>
        inline void Import(const Point3<Q> & b)
        {
            _v[0] = P3ScalarType(b[0]);
            _v[1] = P3ScalarType(b[1]);
            _v[2] = P3ScalarType(b[2]);
        }

        static inline Point3 Construct(const Point3<ScalarType> & b)
        {
            return b;
        }

        inline void SetZero()
        {
            _v[0] = 0;
            _v[1] = 0;
            _v[2] = 0;
        }

        inline const P3ScalarType &X() const
        {
            return _v[0];
        }
        inline const P3ScalarType &Y() const
        {
            return _v[1];
        }
        inline const P3ScalarType &Z() const
        {
            return _v[2];
        }

        inline P3ScalarType &X() { return _v[0]; }
        inline P3ScalarType &Y() { return _v[1]; }
        inline P3ScalarType &Z() { return _v[2]; }

        inline P3ScalarType & operator [] (const int i)
        {
            assert(i >= 0 && i < 3);
            return _v[i];
        }
        inline const P3ScalarType & operator [] (const int i) const
        {
            assert(i >= 0 && i < 3);
            return _v[i];
        }

        inline Point3 operator + (Point3 const &p) const
        {
            return Point3<P3ScalarType>(_v[0] + p._v[0], _v[1] + p._v[1], _v[2] + p._v[2]);
        }
        inline Point3 operator - (Point3 const &p) const
        {
            return Point3<P3ScalarType>(_v[0] - p._v[0], _v[1] - p._v[1], _v[2] - p._v[2]);
        }
        inline Point3 operator * (const P3ScalarType s) const
        {
            return Point3<P3ScalarType>(_v[0] * s, _v[1] * s, _v[2] * s);
        }

        inline Point3 operator / (const P3ScalarType s) const
        {
            return Point3<P3ScalarType>(_v[0] / s, _v[1] / s, _v[2] / s);
        }

        inline Point3 &operator += (Point3 const &p)
        {
            _v[0] += p._v[0];
            _v[1] += p._v[1];
            _v[2] += p._v[2];
            return *this;
        }

        inline Point3 & operator -= (Point3 const & p)
        {
            _v[0] -= p._v[0];
            _v[1] -= p._v[1];
            _v[2] -= p._v[2];
            return *this;
        }

        inline Point3 operator - () const
        {
            return Point3<P3ScalarType>(-_v[0], -_v[1], -_v[2]);
        }

        inline P3ScalarType operator * (Point3 const & p) const
        {
            return (_v[0] * p._v[0] + _v[1] * p._v[1] + _v[2] * p._v[2]);
        }

        inline bool operator == (Point3 const & p) const
        {
            return _v[0] == p._v[0] && _v[1] == p._v[1] && _v[2] == p._v[2];
        }

        inline Point3 & operator *= (const P3ScalarType s)
        {
            _v[0] *= s;
            _v[1] *= s;
            _v[2] *= s;
            return *this;
        }

        inline Point3 & operator /= (const P3ScalarType s)
        {
            _v[0] /= s;
            _v[1] /= s;
            _v[2] /= s;
            return *this;
        }

        inline P3ScalarType dot(const Point3 & p) const { return (*this) * p; }

        inline Point3 operator ^ (Point3 const & p) const
        {
            return Point3 <P3ScalarType>
                (
                    _v[1] * p._v[2] - _v[2] * p._v[1],
                    _v[2] * p._v[0] - _v[0] * p._v[2],
                    _v[0] * p._v[1] - _v[1] * p._v[0]
                    );
        }

        inline P3ScalarType Norm() const
        {
            return math::Sqrt(_v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2]);
        }

        inline P3ScalarType SquaredNorm() const
        {
            return (_v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2]);
        }

        inline Point3 & Normalize()
        {
            P3ScalarType n = P3ScalarType(math::Sqrt(_v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2]));
            if (n > P3ScalarType(0)) { _v[0] /= n; _v[1] /= n; _v[2] /= n; }
            return *this;
        }

    };

    template <class P3ScalarType>
    inline P3ScalarType SquaredNorm(Point3<P3ScalarType> const &p)
    {
        return p.SquaredNorm();
    }

    template <class P3ScalarType>
    inline Point3<P3ScalarType> &Normalize(Point3<P3ScalarType> &p)
    {
        p.Normalize();
        return p;
    }

    template <class P3ScalarType>
    inline P3ScalarType AngleN(Point3<P3ScalarType> const & p1, Point3<P3ScalarType> const & p2)
    {
        P3ScalarType w = p1*p2;
        if (w > 1)
            w = 1;
        else if (w < -1)
            w = -1;
        return (P3ScalarType)acos(w);
    }

    template <class P3ScalarType>
    inline P3ScalarType SquaredDistance(Point3<P3ScalarType> const &p1, Point3<P3ScalarType> const &p2)
    {
        return (p1 - p2).SquaredNorm();
    }


    template <class P3ScalarType>
    inline P3ScalarType Distance(Point3<P3ScalarType> const & p1, Point3<P3ScalarType> const & p2)
    {
        return (p1 - p2).Norm();
    }

    template <class P3ScalarType>
    inline P3ScalarType Norm(Point3<P3ScalarType> const & p)
    {
        return p.Norm();
    }

    template<class P3ScalarType>
    P3ScalarType PSDist(const Point3<P3ScalarType> & p,
        const Point3<P3ScalarType> & v1,
        const Point3<P3ScalarType> & v2,
        Point3<P3ScalarType> & q)
    {
        Point3<P3ScalarType> e = v2 - v1;
        P3ScalarType  t = ((p - v1)*e) / e.SquaredNorm();
        if (t < 0)      t = 0;
        else if (t > 1) t = 1;
        q = v1 + e*t;
        return Distance(p, q);
    }

    typedef Point3<int>	   Point3i;
    typedef Point3<float>  Point3f;    
    
    namespace tri
    {
        template <class UpdateMeshType>
        class UpdateFlags
        {
        public:
            typedef UpdateMeshType MeshType;
            typedef typename MeshType::VertexType     VertexType;
            typedef typename MeshType::VertexIterator VertexIterator;
            typedef typename MeshType::FaceIterator   FaceIterator;

            static void VertexClear(MeshType &m, unsigned int FlagMask = 0xffffffff)
            {
                RequirePerVertexFlags(m);
                int andMask = ~FlagMask;
                for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                    if (!(*vi).IsD()) (*vi).Flags() &= andMask;
            }

            static void VertexSet(MeshType &m, unsigned int FlagMask)
            {
                RequirePerVertexFlags(m);
                for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                    if (!(*vi).IsD()) (*vi).Flags() |= FlagMask;
            }

            static void VertexClearV(MeshType &m) { VertexClear(m, VertexType::VISITED); }
            static void VertexSetV(MeshType &m) { VertexSet(m, VertexType::VISITED); }

            static void FaceBorderFromFF(MeshType &m)
            {
                RequirePerFaceFlags(m);
                RequireFFAdjacency(m);

                for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)if (!(*fi).IsD())
                {
                    for (int j = 0; j < fi->VN(); ++j)
                    {
                        if (face::IsBorder(*fi, j))
                            (*fi).SetB(j);
                        else (*fi).ClearB(j);
                    }
                }
            }

        };

        template <class UpdateMeshType>
        class UpdateTopology
        {
        public:
            typedef UpdateMeshType MeshType;
            typedef typename MeshType::ScalarType     ScalarType;
            typedef typename MeshType::VertexPointer  VertexPointer;
            typedef typename MeshType::VertexIterator VertexIterator;
            typedef typename MeshType::FacePointer    FacePointer;
            typedef typename MeshType::FaceIterator   FaceIterator;

            class PEdge
            {
            public:
                VertexPointer  v[2];
                FacePointer    f;
                int            z;
                bool isBorder;

                PEdge() {}
                PEdge(FacePointer  pf, const int nz)
                {
                    this->Set(pf, nz);
                }
                void Set(FacePointer  pf, const int nz)
                {
                    assert(pf != 0);
                    assert(nz >= 0);
                    assert(nz < pf->VN());

                    v[0] = pf->V(nz);
                    v[1] = pf->V(pf->Next(nz));
                    assert(v[0] != v[1]);

                    if (v[0] > v[1]) std::swap(v[0], v[1]);
                    f = pf;
                    z = nz;
                }

                inline bool operator <  (const PEdge &pe) const
                {
                    if (v[0] < pe.v[0]) return true;
                    else if (v[0] > pe.v[0]) return false;
                    else return v[1] < pe.v[1];
                }

                inline bool operator == (const PEdge &pe) const
                {
                    return v[0] == pe.v[0] && v[1] == pe.v[1];
                }
            };

            static void FillEdgeVector(MeshType &m, std::vector<PEdge> &edgeVec, bool includeFauxEdge = true)
            {
                edgeVec.reserve(m.fn * 3);
                for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
                    if (!(*fi).IsD())
                        for (int j = 0; j < (*fi).VN(); ++j)
                            if (includeFauxEdge || !(*fi).IsF(j))
                                edgeVec.push_back(PEdge(&*fi, j));
            }

            static void FaceFace(MeshType &m)
            {
                RequireFFAdjacency(m);
                if (m.fn == 0) return;

                std::vector<PEdge> e;
                FillEdgeVector(m, e);
                sort(e.begin(), e.end());

                int ne = 0;

                typename std::vector<PEdge>::iterator pe, ps;
                ps = e.begin();
                pe = e.begin();
                do
                {
                    if (pe == e.end() || !(*pe == *ps))
                    {
                        typename std::vector<PEdge>::iterator q, q_next;
                        for (q = ps; q < pe - 1; ++q)
                        {
                            assert((*q).z >= 0);
                            q_next = q;
                            ++q_next;
                            assert((*q_next).z >= 0);
                            assert((*q_next).z < (*q_next).f->VN());
                            (*q).f->FFp(q->z) = (*q_next).f;
                            (*q).f->FFi(q->z) = (*q_next).z;
                        }
                        assert((*q).z >= 0);
                        assert((*q).z < (*q).f->VN());
                        (*q).f->FFp((*q).z) = ps->f;
                        (*q).f->FFi((*q).z) = ps->z;
                        ps = pe;
                        ++ne;
                    }
                    if (pe == e.end()) 
                        break;
                    ++pe;
                }
                while (true);
            }

            static void VertexFace(MeshType &m)
            {
                RequireVFAdjacency(m);

                for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                {
                    (*vi).VFp() = 0;
                    (*vi).VFi() = 0;
                }

                for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
                    if (!(*fi).IsD())
                    {
                        for (int j = 0; j < (*fi).VN(); ++j)
                        {
                            (*fi).VFp(j) = (*fi).V(j)->VFp();
                            (*fi).VFi(j) = (*fi).V(j)->VFi();
                            (*fi).V(j)->VFp() = &(*fi);
                            (*fi).V(j)->VFi() = j;
                        }
                    }
            }

        };

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

    }

    template<class TriangleType>
    typename TriangleType::ScalarType DoubleArea(const TriangleType &t)
    {
        return Norm((t.cP(1) - t.cP(0)) ^ (t.cP(2) - t.cP(0)));
    }

    template <class BoxScalarType>
    class Box3
    {
    public:
        typedef BoxScalarType ScalarType;

        Point3<BoxScalarType> min;
        Point3<BoxScalarType> max;

        inline  Box3() { min.X() = 1; max.X() = -1; min.Y() = 1; max.Y() = -1; min.Z() = 1; max.Z() = -1; }
        inline  Box3(const Box3 & b) { min = b.min; max = b.max; }
        inline  Box3(const Point3<BoxScalarType> & mi, const Point3<BoxScalarType> & ma) { min = mi; max = ma; }
        inline Box3(const Point3<BoxScalarType> & center, const BoxScalarType & radius) {
            min = center - Point3<BoxScalarType>(radius, radius, radius);
            max = center + Point3<BoxScalarType>(radius, radius, radius);
        }
        inline ~Box3() { }

        void Set(const Point3<BoxScalarType> & p)
        {
            min = max = p;
        }

        void SetNull()
        {
            min.X() = 1; max.X() = -1;
            min.Y() = 1; max.Y() = -1;
            min.Z() = 1; max.Z() = -1;
        }

        void Add(Box3<BoxScalarType> const & b)
        {
            if (b.IsNull()) return;
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

        bool IsInEx(Point3<BoxScalarType> const & p) const
        {
            return (
                min.X() <= p.X() && p.X() < max.X() &&
                min.Y() <= p.Y() && p.Y() < max.Y() &&
                min.Z() <= p.Z() && p.Z() < max.Z()
                );
        }

        bool IsNull() const { return min.X() > max.X() || min.Y() > max.Y() || min.Z() > max.Z(); }

        BoxScalarType Diag() const
        {
            return Distance(min, max);
        }

    };

    typedef Box3<int>	 Box3i;
    typedef Box3<float>  Box3f;
    
    template <class SegmentScalarType >
    class Segment3
    {
    public:
        typedef SegmentScalarType ScalarType;
        typedef Point3<SegmentScalarType> PointType;
        typedef Segment3<SegmentScalarType> SegmentType;

    private:
        PointType _p0, _p1;

    public:
        inline const PointType &P0() const { return _p0; }
        inline const PointType &P1() const { return _p1; }

        Segment3() {};
        Segment3(const PointType &a, const PointType &b) { _p0 = a; _p1 = b; };

        inline PointType MidPoint() const
        {
            return (_p0 + _p1) / ScalarType(2.0);
        }

    };

    typedef Segment3<float>  Segment3f;

    template <class ScalarType>
    void SegmentPointDistance(Segment3<ScalarType> s,
        const Point3<ScalarType> & p,
        Point3< ScalarType > &clos,
        ScalarType &dist)
    {
        SegmentPointSquaredDistance(s, p, clos, dist);
        dist = sqrt(dist);
    }

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


    template <class T, bool NORM = true> class Plane3 {
    public:
        typedef T ScalarType;
        typedef Point3<T> PointType;

    private:
        ScalarType _offset;
        PointType _dir;

    public:
        Plane3() {}
        Plane3(const ScalarType &dist, const PointType &dir) { Set(dist, dir); }

        const ScalarType &Offset() const { return _offset; }
        ScalarType &Offset() { return _offset; }
        void SetOffset(const ScalarType &o) { _offset = o; }

        const PointType &Direction() const { return _dir; }

        void Normalize() {
            _dir.Normalize();
        }

        inline void Init(const PointType &p0, const PointType &norm) {
            _dir = norm;
            if (NORM) Normalize();
            _offset = p0.dot(_dir);
        }
    };

    typedef Plane3<float>  Plane3f;

    template<class T> T SignedDistancePlanePoint(const Plane3<T, true> & plane, const Point3<T> & point)
    {
        return plane.Direction().dot(point) - plane.Offset();
    }

    template <class SCALARTYPE>
    class BasicGrid
    {
    public:
        typedef SCALARTYPE ScalarType;
        typedef Box3<ScalarType> Box3x;
        typedef Point3<ScalarType> CoordType;
        typedef BasicGrid<SCALARTYPE> GridType;

        Box3x bbox;

        CoordType dim;
        Point3i siz;
        CoordType voxel;

        inline void PToIP(const CoordType & p, Point3i &pi) const
        {
            CoordType t = p - bbox.min;
            pi[0] = int(t[0] / voxel[0]);
            pi[1] = int(t[1] / voxel[1]);
            pi[2] = int(t[2] / voxel[2]);
        }

        inline void BoxToIBox(const Box3x & b, Box3i & ib) const
        {
            PToIP(b.min, ib.min);
            PToIP(b.max, ib.max);
        }
    };

    template<class scalar_type>
    void BestDim(const __int64 elems, const Point3<scalar_type> & size, Point3i & dim)
    {
        const __int64 mincells = 1;
        const double GFactor = 1;
        double diag = size.Norm();
        double eps = diag*1e-4;

        assert(elems>0);
        assert(size[0] >= 0.0);
        assert(size[1] >= 0.0);
        assert(size[2] >= 0.0);


        __int64 ncell = (__int64)(elems*GFactor);
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

    template <class OBJTYPE, class SCALARTYPE>
    class SpatialIndex {
    public:
        typedef SpatialIndex<OBJTYPE, SCALARTYPE> ClassType;
        typedef OBJTYPE ObjType;
        typedef SCALARTYPE ScalarType;
        typedef ObjType * ObjPtr;
        typedef Point3<ScalarType> CoordType;
        typedef Box3<ScalarType> BoxType;

        template <class OBJMARKER, class OBJPTRCONTAINER>
        unsigned int GetInBox(OBJMARKER & _marker, const BoxType _bbox, OBJPTRCONTAINER & _objectPtrs) {
            assert(0);
            (void)_marker;
            (void)_bbox;
            (void)_objectPtrs;
            return (0);
        }
    };
    
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

        class Link
        {
        public:
            inline Link() {};
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
            ObjPtr t;
            int i;
        };

        typedef Link* Cell;
        typedef Cell CellIterator;

        std::vector<Link>   links;

        std::vector<Cell> grid;


        bool Empty() const { return links.empty(); }

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

            ScalarType infl = _bbox.Diag() / _size;
            _bbox.min -= tessellation::Point3<FLT>(infl, infl, infl);
            _bbox.max += tessellation::Point3<FLT>(infl, infl, infl);

            Set(_oBegin, _oEnd, _bbox, _size);
        }

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

        template <class OBJITER>
        inline void Set(const OBJITER & _oBegin, const OBJITER & _oEnd, const Box3x &_bbox, Point3i _siz)
        {
            OBJITER i;

            this->bbox = _bbox;
            this->siz = _siz;
            
            this->dim = this->bbox.max - this->bbox.min;
            this->voxel[0] = this->dim[0] / this->siz[0];
            this->voxel[1] = this->dim[1] / this->siz[1];
            this->voxel[2] = this->dim[2] / this->siz[2];

            grid.resize(this->siz[0] * this->siz[1] * this->siz[2] + 1);

            links.clear();
            for (i = _oBegin; i != _oEnd; ++i)
            {
                Box3x bb;
                (*i).GetBBox(bb);
                bb.Intersect(this->bbox);
                if (!bb.IsNull())
                {
                    Box3i ib;
                    this->BoxToIBox(bb, ib);
                    int x, y, z;
                    for (z = ib.min[2]; z <= ib.max[2]; ++z)
                    {
                        int bz = z*this->siz[1];
                        for (y = ib.min[1]; y <= ib.max[1]; ++y)
                        {
                            int by = (y + bz)*this->siz[0];
                            for (x = ib.min[0]; x <= ib.max[0]; ++x)
                                links.push_back(Link(&(*i), by + x));
                        }
                    }
                }
            }

            links.push_back(Link(NULL, int(grid.size()) - 1));

            sort(links.begin(), links.end());

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
        ObjPtr GetClosest(OBJPOINTDISTFUNCTOR & _getPointDistance, OBJMARKER & _marker,
            const typename OBJPOINTDISTFUNCTOR::QueryType & _p, const ScalarType & _maxDist, ScalarType & _minDist, CoordType & _closestPt)
        {
            return (tessellation::GridClosest<GridPtrType, OBJPOINTDISTFUNCTOR, OBJMARKER>(*this, _getPointDistance, _marker, _p, _maxDist, _minDist, _closestPt));
        }

    };


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
                    if (_getPointDistance((**l), _p_obj, _minDist, t_res))
                    {
                        winner = elem;
                        _closestPt = t_res;
                        newradius = _minDist;
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

            Si.BoxToIBox(boxtodo, iboxtodo);
            iboxtodo.Intersect(ibox);
            if (!boxtodo.IsNull())
            {
                for (ix = iboxtodo.min[0]; ix <= iboxtodo.max[0]; ix++)
                    for (iy = iboxtodo.min[1]; iy <= iboxtodo.max[1]; iy++)
                        for (iz = iboxtodo.min[2]; iz <= iboxtodo.max[2]; iz++)
                            if (ix<iboxdone.min[0] || ix>iboxdone.max[0] ||
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
}
