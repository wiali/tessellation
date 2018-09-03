#pragma  once

#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace tessellation
{
    namespace math
    {
        const float EPSILON = 0.001f;

        template <class SCALAR>
        inline SCALAR  Clamp(const SCALAR &val, const SCALAR &minval, const SCALAR &maxval)
        {
            if (val < minval) return minval;
            if (val > maxval) return maxval;
            return val;
        }

        inline float Sqrt(const float v) { return sqrtf(v); }

        inline bool isZero(const float v)
        {
            return fabs(v) < EPSILON;
        }

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
                else return c; // if c<a then c is smaller than b...
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


    // chain with 2 template arguments on the derivers
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
        /// The only data member. Hidden to user.
        P2ScalarType _v[2];
    public:
        /// the scalar type
        typedef P2ScalarType ScalarType;
        enum { Dimension = 2 };

        /// empty constructor (does nothing)
        inline Point2() { }
        /// x,y constructor
        inline Point2(const ScalarType nx, const ScalarType ny)
        {
            _v[0] = nx;
            _v[1] = ny;
        }

        inline const ScalarType &X() const { return _v[0]; }
        inline const ScalarType &Y() const { return _v[1]; }
        inline ScalarType &X() { return _v[0]; }
        inline ScalarType &Y() { return _v[1]; }

        //@{
        /** @name Linearity for 2d points (operators +, -, *, /, *= ...) **/
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

        /// normalizes, and returns itself as result
        inline Point2 & Normalize(void)
        {
            ScalarType n = math::Sqrt(_v[0] * _v[0] + _v[1] * _v[1]);
            if (n > 0.0) { _v[0] /= n;	_v[1] /= n; }
            return *this;
        }

    }; // end class definition

    typedef Point2<float>  Point2f;


    template <class T> class Box3;

    template <class P3ScalarType> class Point3
    {
    protected:
        /// The only data member. Hidden to user.
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

        //@}
        //@{

        /** @name Classical overloading of operators
        Note
        **/

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
        /// Dot product


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
        /// Cross product
        inline Point3 operator ^ (Point3 const & p) const
        {
            return Point3 <P3ScalarType>
                (
                    _v[1] * p._v[2] - _v[2] * p._v[1],
                    _v[2] * p._v[0] - _v[0] * p._v[2],
                    _v[0] * p._v[1] - _v[1] * p._v[0]
                    );
        }

        // Norme
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

    }; // end class definition

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

    /// Point(p) Edge(v1-v2) dist, q is the point in v1-v2 with min dist
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

    template <class T> class Point4
    {
    public:
        /// The only data member. Hidden to user.
        T _v[4];

    public:
        typedef T ScalarType;
        enum { Dimension = 4 };

        //@{

        /** @name Standard Constructors and Initializers
         No casting operators have been introduced to avoid automatic unattended (and costly) conversion between different point types
         **/

        inline Point4() { }
        inline Point4(const T nx, const T ny, const T nz, const T nw)
        {
            _v[0] = nx;
            _v[1] = ny;
            _v[2] = nz;
            _v[3] = nw;
        }


        /** @name Data Access.
         access to data is done by overloading of [] or explicit naming of coords (x,y,z,w)
          **/
        inline const T &operator [] (const int i) const
        {
            assert(i >= 0 && i < 4);
            return _v[i];
        }
        inline T &operator [] (const int i)
        {
            assert(i >= 0 && i < 4);
            return _v[i];
        }

        inline T const *V() const
        {
            return _v;
        }
        inline T *V()
        {
            return _v;
        }

    }; // end class definition


    template <class T>
    class Color4 : public Point4<T>
    {
        typedef Point4<T> Base;
    public:
        /// Constant for storing standard colors.
        /// Each color is stored in a simple in so that the bit pattern match with the one of Color4b.
        enum ColorConstant
        {
            Black = 0xff000000,
            Gray = 0xff808080,
            White = 0xffffffff,

            Red = 0xff0000ff,
            Green = 0xff00ff00,
            Blue = 0xffff0000,

            Cyan = 0xffffff00,
            Yellow = 0xff00ffff,
            Magenta = 0xffff00ff,

            LightGray = 0xffc0c0c0,
            LightRed = 0xff8080ff,
            LightGreen = 0xff80ff80,
            LightBlue = 0xffff8080,

            DarkGray = 0xff404040,
            DarkRed = 0xff000040,
            DarkGreen = 0xff004000,
            DarkBlue = 0xff400000
        };

        inline Color4(const T nx, const T ny, const T nz, const T nw) : Point4<T>(nx, ny, nz, nw) {}
        inline Color4(const Point4<T> &c) : Point4<T>(c) {}
        inline Color4() {}
        inline Color4(ColorConstant cc);
        inline Color4(unsigned int cc);


        template <class ScalarInterpType>
        inline void lerp(const Color4 &c0, const Color4 &c1, const ScalarInterpType x)
        {
            assert(x >= 0);
            assert(x <= 1);

            (*this)[0] = (T)(c1.V()[0] * x + c0.V()[0] * (1.0f - x));
            (*this)[1] = (T)(c1.V()[1] * x + c0.V()[1] * (1.0f - x));
            (*this)[2] = (T)(c1.V()[2] * x + c0.V()[2] * (1.0f - x));
            (*this)[3] = (T)(c1.V()[3] * x + c0.V()[3] * (1.0f - x));
        }
    }; /// END CLASS ///////////////////

    template<>
    inline Color4<unsigned char>::Color4(Color4<unsigned char>::ColorConstant cc)
    {
        *((int *)this) = cc;
    }


    typedef Color4<unsigned char>  Color4b;



    namespace tri
    {
        /// \ingroup trimesh

        /// \headerfile flag.h vcg/complex/algorithms/update/flag.h

        /// \brief Management, updating and computation of per-vertex and per-face flags (like border flags).

        /**
        This class is used to compute or update some of the flags that can be stored in the mesh components. For now just Border flags (e.g. the flag that tells if a given edge of a face belong to a border of the mesh or not).
        */

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
            /**
            \warning Obviously it assumes that the topology has been correctly computed (see: UpdateTopology::FaceFace )
            */
            static void FaceBorderFromFF(MeshType &m)
            {
                RequirePerFaceFlags(m);
                RequireFFAdjacency(m);

                for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)if (!(*fi).IsD())
                        for (int j = 0; j < fi->VN(); ++j)
                        {
                            if (face::IsBorder(*fi, j)) (*fi).SetB(j);
                            else (*fi).ClearB(j);
                        }
            }

        }; // end class

        /// \brief Generation of per-vertex and per-face topological information.

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


            /// \headerfile topology.h vcg/complex/algorithms/update/topology.h

            /// \brief Auxiliairy data structure for computing face face adjacency information.
            /**
            It identifies and edge storing two vertex pointer and a face pointer where it belong.
            */

            class PEdge
            {
            public:

                VertexPointer  v[2];  // the two Vertex pointer are ordered!
                FacePointer    f;     // the face where this edge belong
                int            z;     // index in [0..2] of the edge of the face
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
                    assert(v[0] != v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

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
                /// Convert from edge barycentric coord to the face baricentric coord a point on the current edge.
                /// Face barycentric coordinates are relative to the edge face.
                inline Point3<ScalarType> EdgeBarycentricToFaceBarycentric(ScalarType u) const
                {
                    Point3<ScalarType> interp(0, 0, 0);
                    interp[this->z] = u;
                    interp[(this->z + 1) % 3] = 1.0f - u;
                    return interp;
                }
            };

            /// Fill a vector with all the edges of the mesh.
            /// each edge is stored in the vector the number of times that it appears in the mesh, with the referring face.
            /// optionally it can skip the faux edges (to retrieve only the real edges of a triangulated polygonal mesh)

            static void FillEdgeVector(MeshType &m, std::vector<PEdge> &edgeVec, bool includeFauxEdge = true)
            {
                edgeVec.reserve(m.fn * 3);
                for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
                    if (!(*fi).IsD())
                        for (int j = 0; j < (*fi).VN(); ++j)
                            if (includeFauxEdge || !(*fi).IsF(j))
                                edgeVec.push_back(PEdge(&*fi, j));
            }

            /// \brief Update the Face-Face topological relation by allowing to retrieve for each face what other faces shares their edges.
            static void FaceFace(MeshType &m)
            {
                RequireFFAdjacency(m);
                if (m.fn == 0) return;

                std::vector<PEdge> e;
                FillEdgeVector(m, e);
                sort(e.begin(), e.end());							// Lo ordino per vertici

                int ne = 0;											// Numero di edge reali

                typename std::vector<PEdge>::iterator pe, ps;
                ps = e.begin();
                pe = e.begin();
                //for(ps = e.begin(),pe=e.begin();pe<=e.end();++pe)	// Scansione vettore ausiliario
                do
                {
                    if (pe == e.end() || !(*pe == *ps))					// Trovo blocco di edge uguali
                    {
                        typename std::vector<PEdge>::iterator q, q_next;
                        for (q = ps; q < pe - 1; ++q)						// Scansione facce associate
                        {
                            assert((*q).z >= 0);
                            //assert((*q).z< 3);
                            q_next = q;
                            ++q_next;
                            assert((*q_next).z >= 0);
                            assert((*q_next).z < (*q_next).f->VN());
                            (*q).f->FFp(q->z) = (*q_next).f;				// Collegamento in lista delle facce
                            (*q).f->FFi(q->z) = (*q_next).z;
                        }
                        assert((*q).z >= 0);
                        assert((*q).z < (*q).f->VN());
                        (*q).f->FFp((*q).z) = ps->f;
                        (*q).f->FFi((*q).z) = ps->z;
                        ps = pe;
                        ++ne;										// Aggiorno il numero di edge
                    }
                    if (pe == e.end()) 
                        break;
                    ++pe;
                }
                while (true);
            }

            /// \brief Update the Vertex-Face topological relation.
            /**
            The function allows to retrieve for each vertex the list of faces sharing this vertex.
            After this call all the VF component are initialized. Isolated vertices have a null list of faces.
            \sa vcg::vertex::VFAdj
            \sa vcg::face::VFAdj
            */

            static void VertexFace(MeshType &m)
            {
                RequireVFAdjacency(m);

                for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                {
                    (*vi).VFp() = 0;
                    (*vi).VFi() = 0; // note that (0,-1) means uninitiazlied while 0,0 is the valid initialized values for isolated vertices.
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

        }; // end class

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

    } // end namespace

    template<class TriangleType>
    typename TriangleType::ScalarType DoubleArea(const TriangleType &t)
    {
        return Norm((t.cP(1) - t.cP(0)) ^ (t.cP(2) - t.cP(0)));
    }


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

        bool IsNull() const { return min.X() > max.X() || min.Y() > max.Y() || min.Z() > max.Z(); }

        BoxScalarType Diag() const
        {
            return Distance(min, max);
        }

    }; // end class definition


    typedef Box3<int>	 Box3i;
    typedef Box3<float>  Box3f;
    

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

        /// The empty constructor
        Segment3() {};
        /// The (a,b) constructor
        Segment3(const PointType &a, const PointType &b) { _p0 = a; _p1 = b; };

        /// return the middle point
        inline PointType MidPoint() const
        {
            return (_p0 + _p1) / ScalarType(2.0);
        }

    }; // end class definition

    typedef Segment3<float>  Segment3f;


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


        /// Function to normalize direction
        void Normalize() {
            _dir.Normalize();
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
}
