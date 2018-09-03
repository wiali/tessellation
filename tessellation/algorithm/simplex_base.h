#pragma once

namespace tessellation
{
    struct AllTypes
    {
        struct AVertexType {};
        struct AFaceType {};

    };

    template<class T = float, int NMAX = 1>
    class TexCoord2
    {
    public:
        typedef Point2<T>  PointType;
        typedef T ScalarType;


    private:
        PointType _t[NMAX];
        short     _n[NMAX];
    public:

        TexCoord2(T u, T v)
        {
            if (NMAX > 0) _n[0] = 0;
            _t[0][0] = u;
            _t[0][1] = v;
        }
        TexCoord2() {  }

        /**/inline short      &n()
        {
            return _n[0];
        }
        /**/inline short n() const
        {
            return _n[0];
        }

        /**/inline Point2<T> &t()
        {
            return _t[0];
        }
        /**/inline Point2<T> t() const
        {
            return _t[0];
        }

        inline const PointType &P() const { return _t[0]; }
        inline PointType &P() { return _t[0]; }

        inline const PointType &P(const int i) const { assert(i > 0 && i < NMAX); return _t[i]; }
        inline PointType &P(const int i) { assert(i > 0 && i < NMAX); return _t[i]; }

        /**/inline T & u() { return _t[0][0]; }
        /**/inline T & v() { return _t[0][1]; }
        /**/inline const T & u() const { return _t[0][0]; }
        /**/inline const T & v() const { return _t[0][1]; }

        enum { n_coords = NMAX };
    };

    typedef TexCoord2<float>  TexCoord2f;

    namespace vertex
    {
        /** \addtogroup VertexComponentGroup
        @{

        /*------------------------- EMPTY CORE COMPONENTS -----------------------------------------*/

        template <class TT> class EmptyCore : public TT
        {
        public:
            typedef int FlagType;
            int &Flags()
            {
                assert(0);
                static int dummyflags(0);
                return dummyflags;
            }
            int cFlags() const
            {
                return 0;
            }
            static bool HasFlags()
            {
                return false;
            }

            typedef tessellation::Point3f CoordType;
            typedef CoordType::ScalarType      ScalarType;
            CoordType &P()
            {
                assert(0);
                static CoordType coord(0, 0, 0);
                return coord;
            }
            CoordType cP() const
            {
                assert(0);
                static CoordType coord(0, 0, 0);
                assert(0);
                return coord;
            }
            static bool HasCoord()
            {
                return false;
            }
            inline bool IsCoordEnabled() const
            {
                return TT::VertexType::HasCoord();
            }

            typedef tessellation::Point3f NormalType;
            NormalType &N()
            {
                assert(0);
                static NormalType dummy_normal(0, 0, 0);
                return dummy_normal;
            }
            NormalType cN() const
            {
                assert(0);
                static NormalType dummy_normal(0, 0, 0);
                return dummy_normal;
            }
            static bool HasNormal()
            {
                return false;
            }
            inline bool IsNormalEnabled() const
            {
                return TT::VertexType::HasNormal();
            }

            typedef tessellation::TexCoord2<float, 1> TexCoordType;
            TexCoordType &T()
            {
                static TexCoordType dummy_texcoord;
                assert(0);
                return dummy_texcoord;
            }
            TexCoordType cT() const
            {
                static TexCoordType dummy_texcoord;
                assert(0);
                return dummy_texcoord;
            }
            static bool HasTexCoord()
            {
                return false;
            }
            inline bool IsTexCoordEnabled() const
            {
                return TT::VertexType::HasTexCoord();
            }

            typename TT::FacePointer &VFp()
            {
                static typename TT::FacePointer fp = 0;
                assert(0);
                return fp;
            }
            typename TT::FacePointer cVFp() const
            {
                static typename TT::FacePointer fp = 0;
                assert(0);
                return fp;
            }
            int &VFi()
            {
                static int z = -1;
                assert(0);
                return z;
            }
            int cVFi() const
            {
                static int z = -1;
                assert(0);
                return z;
            }


            int &VEi() { static int z = -1; return z; }
            int cVEi() const { static int z = -1; return z; }

            bool IsNull() const
            {
                return true;
            }
            static bool HasVFAdjacency()
            {
                return false;
            }

            bool IsVFInitialized() const { return static_cast<const typename TT::VertexType *>(this)->cVFi() != -1; }
            void VFClear() {
                if (IsVFInitialized()) {
                    static_cast<typename TT::VertexPointer>(this)->VFp() = 0;
                    static_cast<typename TT::VertexPointer>(this)->VFi() = -1;
                }
            }

            bool IsVEInitialized() const { return static_cast<const typename TT::VertexType *>(this)->cVEi() != -1; }
            void VEClear() {
                if (IsVEInitialized()) {
                    static_cast<typename TT::VertexPointer>(this)->VEi() = -1;
                }
            }

            template < class RightValueType>
            void ImportData(const RightValueType  & /*rVert*/) {
                //			TT::ImportData( rVert);
            }

            static bool HasVEAdjacency() { return false; }
        };

        /*-------------------------- COORD ----------------------------------------*/
        /*! \brief \em Generic Component: \b Geometric \b Position of the vertex
        Templated on the coordinate class. In practice you use one of the two specialized class Coord3f and Coord3d
        You can access to the coordinate of a vertex by mean of the P(),cP() member functions.
        */
        template <class A, class T> class Coord : public T
        {
        public:
            typedef A CoordType;
            typedef typename A::ScalarType      ScalarType;
            /// Return a const reference to the coordinate of the vertex
            inline const CoordType &P() const
            {
                return _coord;
            }
            /// Return a reference to the coordinate of the vertex
            inline       CoordType &P()
            {
                return _coord;
            }
            /// Return a const reference to the coordinate of the vertex
            inline       CoordType cP() const
            {
                return _coord;
            }

            template < class RightValueType>
            void ImportData(const RightValueType   &rVert)
            {
                if (rVert.IsCoordEnabled()) P().Import(rVert.cP());
                T::ImportData(rVert);
            }
            static bool HasCoord()
            {
                return true;
            }

        private:
            CoordType _coord;
        };

        /// Specialized Coord Component in floating point precision.
        template <class T> class Coord3f : public Coord<Point3f, T> {
        public:	static void Name(std::vector<std::string> & name) { name.push_back(std::string("Coord3f")); T::Name(name); }
        };

        /*-------------------------- NORMAL ----------------------------------------*/
        /*! \brief \em Generic Component: \b %Normal of the vertex

        Templated on the Point3 class used to store the normal.
        In practice you use one of the two specialized class Normal3f and Normal3d.

        You can access to the normal of a vertex by mean of the N(),cN() member functions.

        \note Many algorithms assume that, for sake of precision coherence,
        the type of the normal is the same with respect to the type coord component.
        */

        template <class A, class T> class Normal : public T
        {
        public:
            typedef A NormalType;
            /// Return a const reference to the normal of the vertex
            inline const NormalType &N() const
            {
                return _norm;
            }
            /// Return a  reference to the normal of the vertex
            inline       NormalType &N()
            {
                return _norm;
            }
            /// Return a const reference to the normal of the vertex
            inline       NormalType cN() const
            {
                return _norm;
            }
            template < class RightValueType>
            void ImportData(const RightValueType   &rVert)
            {
                if (rVert.IsNormalEnabled())  N().Import(rVert.cN());
                T::ImportData(rVert);
            }
            static bool HasNormal()
            {
                return true;
            }

        private:
            NormalType _norm;
        };

        /// Specialized Normal component in floating point precision.
        template <class T> class Normal3f : public Normal<Point3f, T> {
        public:	static void Name(std::vector<std::string> & name) { name.push_back(std::string("Normal3f")); T::Name(name); }
        };

        template <class A, class TT> class TexCoord : public TT
        {
        public:
            typedef A TexCoordType;

            /// Return a const reference to the Texture Coordinate
            const TexCoordType &T() const
            {
                return _t;
            }
            TexCoordType &T()
            {
                return _t;
            }
            TexCoordType cT() const
            {
                return _t;
            }
            template < class RightValueType>
            void ImportData(const RightValueType   &rVert)
            {
                if (rVert.IsTexCoordEnabled())  T() = rVert.cT();
                TT::ImportData(rVert);
            }
            static bool HasTexCoord()
            {
                return true;
            }
            static void Name(std::vector<std::string> &name)
            {
                name.push_back(std::string("TexCoord"));
                TT::Name(name);
            }

        private:
            TexCoordType _t;
        };

        /// Specialized Texture component in floating point precision.
        template <class TT> class TexCoord2f : public TexCoord<TexCoord2<float, 1>, TT>
        {
        public:
            static void Name(std::vector<std::string> &name)
            {
                name.push_back(std::string("TexCoord2f"));
                TT::Name(name);
            }
        };
        /*------------------------- FLAGS -----------------------------------------*/
        /*! \brief \em Component: Per vertex \b Flags

        This component stores a 32 bit array of bit flags.
        These bit flags are used for keeping track of selection, deletion, visiting etc. \sa \ref flags for more details on common uses of flags.
        */

        template <class T> class BitFlags : public T
        {
        public:
            BitFlags()
            {
                _flags = 0;
            }
            typedef int FlagType;
            inline const int &Flags() const
            {
                return _flags;
            }
            inline       int &Flags()
            {
                return _flags;
            }
            inline       int cFlags() const
            {
                return _flags;
            }

            static bool HasFlags() { return true; }

        private:
            int  _flags;
        };

        template <class VALUE_TYPE>
        class vector_ocf : public std::vector<VALUE_TYPE> {
            typedef std::vector<VALUE_TYPE> BaseType;
            typedef typename vector_ocf<VALUE_TYPE>::iterator ThisTypeIterator;

        public:
            vector_ocf() :std::vector<VALUE_TYPE>()
            {
                MarkEnabled = false;
                NormalEnabled = false;
                TexCoordEnabled = false;
                VFAdjacencyEnabled = false;
            }

            ////////////////////////////////////////
            // All the standard methods of std::vector that can change the reallocation are
            // redefined in order to manage the additional data.
            void push_back(const VALUE_TYPE & v)
            {
                BaseType::push_back(v);
                BaseType::back()._ovp = this;
                if (MarkEnabled)          MV.push_back(0);
                if (NormalEnabled)        NV.push_back(typename VALUE_TYPE::NormalType());
                if (TexCoordEnabled)      TV.push_back(typename VALUE_TYPE::TexCoordType());
                if (VFAdjacencyEnabled)   AV.push_back(VFAdjType());
            }

            void pop_back();

            virtual void resize(size_t _size)
            {
                const size_t oldsize = BaseType::size();
                BaseType::resize(_size);
                if (oldsize < _size) {
                    ThisTypeIterator firstnew = BaseType::begin();
                    advance(firstnew, oldsize);
                    _updateOVP(firstnew, (*this).end());
                }
                if (MarkEnabled)          MV.resize(_size);
                if (NormalEnabled)        NV.resize(_size);
                if (TexCoordEnabled)      TV.resize(_size);
                if (VFAdjacencyEnabled)   AV.resize(_size, VFAdjType::Zero());
            }

            void reserve(size_t _size)
            {
                BaseType::reserve(_size);
                if (MarkEnabled)         MV.reserve(_size);
                if (NormalEnabled)       NV.reserve(_size);
                if (TexCoordEnabled)     TV.reserve(_size);
                if (VFAdjacencyEnabled)  AV.reserve(_size);
            }

            void _updateOVP(ThisTypeIterator lbegin, ThisTypeIterator lend)
            {
                ThisTypeIterator vi;
                for (vi = lbegin; vi != lend; ++vi)
                    (*vi)._ovp = this;
            }

            bool IsMarkEnabled() const { return MarkEnabled; }
            void EnableMark() {
                assert(VALUE_TYPE::HasMarkOcf());
                MarkEnabled = true;
                MV.resize((*this).size(), 0);
            }
            void DisableMark() {
                assert(VALUE_TYPE::HasMarkOcf());
                MarkEnabled = false;
                MV.clear();
            }

            bool IsNormalEnabled() const { return NormalEnabled; }
            void EnableNormal() {
                assert(VALUE_TYPE::HasNormalOcf());
                NormalEnabled = true;
                NV.resize((*this).size());
            }
            void DisableNormal() {
                assert(VALUE_TYPE::HasNormalOcf());
                NormalEnabled = false;
                NV.clear();
            }

            bool IsVFAdjacencyEnabled() const { return VFAdjacencyEnabled; }
            void EnableVFAdjacency() {
                assert(VALUE_TYPE::HasVFAdjacencyOcf());
                VFAdjacencyEnabled = true;
                AV.resize((*this).size(), VFAdjType::Zero());
            }
            void DisableVFAdjacency() {
                assert(VALUE_TYPE::HasVFAdjacencyOcf());
                VFAdjacencyEnabled = false;
                AV.clear();
            }

            bool IsTexCoordEnabled() const { return TexCoordEnabled; }
            void EnableTexCoord() {
                assert(VALUE_TYPE::HasTexCoordOcf());
                TexCoordEnabled = true;
                TV.resize((*this).size());
            }
            void DisableTexCoord() {
                assert(VALUE_TYPE::HasTexCoordOcf());
                TexCoordEnabled = false;
                TV.clear();
            }

            struct VFAdjType {
                VFAdjType() :_fp(0), _zp(-1) {}
                VFAdjType(typename VALUE_TYPE::FacePointer fp, int zp) :_fp(fp), _zp(zp) {}
                typename VALUE_TYPE::FacePointer _fp;
                int _zp;
                static VFAdjType Zero() { return VFAdjType(0, -1); }
                bool IsNull() const { return (_zp == -1); }
            };

        public:
            std::vector<int> MV;
            std::vector<typename VALUE_TYPE::NormalType> NV;
            std::vector<typename VALUE_TYPE::TexCoordType> TV;
            std::vector<struct VFAdjType> AV;

            bool MarkEnabled;
            bool NormalEnabled;
            bool TexCoordEnabled;
            bool VFAdjacencyEnabled;
        };

        ///*-------------------------- TEXTURE  ----------------------------------*/

        template <class A, class TT> class TexCoordOcf : public TT {
        public:
            typedef A TexCoordType;
            const TexCoordType &T() const { assert((*this).Base().TexCoordEnabled); return (*this).Base().TV[(*this).Index()]; }
            TexCoordType &T() { assert((*this).Base().TexCoordEnabled); return (*this).Base().TV[(*this).Index()]; }
            TexCoordType cT() const { assert((*this).Base().TexCoordEnabled); return (*this).Base().TV[(*this).Index()]; }
            template < class RightVertexType>
            void ImportData(const RightVertexType & rightV)
            {
                if ((*this).IsTexCoordEnabled() && rightV.IsTexCoordEnabled()) // copy the data only if they are enabled in both vertices
                    T() = rightV.cT();
                TT::ImportData(rightV);
            }
            inline bool IsTexCoordEnabled()     const { return this->Base().IsTexCoordEnabled(); }
            static bool HasTexCoord() { return true; }
            static bool HasTexCoordOcf() { assert(!TT::HasTexCoordOcf()); return true; }
        };

        template <class T> class TexCoordfOcf : public TexCoordOcf<TexCoord2<float, 1>, T> {
        public: static void Name(std::vector<std::string> & name) { name.push_back(std::string("TexCoordfOcf")); T::Name(name); }
        };

        ///*-------------------------- InfoOpt  ----------------------------------*/

        template < class T> class InfoOcf : public T {
        public:
            // You should never ever try to copy a vertex that has OCF stuff.
            // use ImportData function.
            inline InfoOcf &operator=(const InfoOcf & /*other*/) {
                assert(0); return *this;
            }

            vector_ocf<typename T::VertexType> &Base() const { return *_ovp; }

            inline int Index() const {
                typename  T::VertexType const *tp = static_cast<typename T::VertexType const*>(this);
                int tt2 = tp - &*(_ovp->begin());
                return tt2;
            }
        public:
            vector_ocf<typename T::VertexType> *_ovp;

            static bool HasNormalOcf() { return false; }
            static bool HasMarkOcf() { return false; }
            static bool HasTexCoordOcf() { return false; }
            static bool HasVFAdjacencyOcf() { return false; }
        };

        /*----------------------------- VFADJ ------------------------------*/


        template <class T> class VFAdjOcf : public T {
        public:
            typename T::FacePointer &VFp() {
                assert((*this).Base().VFAdjacencyEnabled);
                return (*this).Base().AV[(*this).Index()]._fp;
            }
            typename T::FacePointer cVFp() const {
                if (!(*this).Base().VFAdjacencyEnabled) return 0;
                else return (*this).Base().AV[(*this).Index()]._fp;
            }

            int &VFi() {
                assert((*this).Base().VFAdjacencyEnabled);
                return (*this).Base().AV[(*this).Index()]._zp;
            }
            int cVFi() const {
                if (!(*this).Base().VFAdjacencyEnabled) return -1;
                return (*this).Base().AV[(*this).Index()]._zp;
            }
            template <class RightVertexType>
            void ImportData(const RightVertexType & rightV)
            {
                T::ImportData(rightV);
            }

            static bool HasVFAdjacency() { return true; }
            static bool HasVFAdjacencyOcf() { return true; }
            bool IsVFAdjacencyEnabled(const typename T::VertexType *vp) { return vp->Base().VFAdjacencyEnabled; }

            static void Name(std::vector<std::string> & name) { name.push_back(std::string("VFAdjOcf")); T::Name(name); }
        private:
        };

        /** @} */   // End Doxygen VertexComponentGroup
    } // end namespace vert

    /* The base class form which we start to add our components.
    it has the empty definition for all the standard members (coords, color flags)
    Note:
    in order to avoid both virtual classes and ambiguous definitions all
    the subsequent overrides must be done in a sequence of derivation.

    In other words we cannot derive and add in a single derivation step
    (with multiple ancestor), both the real (non-empty) normal and color but
    we have to build the type a step a time (deriving from a single ancestor at a time).

     The Real Big Vertex class;

    The class __VertexArityMax__ is the one that is the Last to be derived,
    and therefore is the only one to know the real members
    (after the many overrides) so all the functions with common behaviour
    using the members defined in the various Empty/nonEmpty component classes
    MUST be defined here.

    I.e. IsD() that uses the overridden Flags() member must be defined here.

    */

    template < class UserTypes,
        template <typename> class A, template <typename> class B,
        template <typename> class C, template <typename> class D,
        template <typename> class E, template <typename> class F>
    class VertexArityMax : public Arity6<vertex::EmptyCore<UserTypes>, A, B, C, D, E, F>
    {

        // ----- Flags stuff -----
    public:
        enum
        {

            DELETED = 0x0001,		// This bit indicate that the vertex is deleted from the mesh
            NOTREAD = 0x0002,		// This bit indicate that the vertex of the mesh is not readable
            NOTWRITE = 0x0004,		// This bit indicate that the vertex is not modifiable
            MODIFIED = 0x0008,		// This bit indicate that the vertex is modified
            VISITED = 0x0010,		// This bit can be used to mark the visited vertex
            SELECTED = 0x0020,		// This bit can be used to select
            BORDER = 0x0100,    // Border Flag
            USER0 = 0x0200			// First user bit
        };

        bool IsD() const
        {
            return (this->cFlags() & DELETED) != 0;    ///  checks if the vertex is deleted
        }
        bool IsR() const
        {
            return (this->cFlags() & NOTREAD) == 0;    ///  checks if the vertex is readable
        }
        bool IsW() const
        {
            return (this->cFlags() & NOTWRITE) == 0;    ///  checks if the vertex is modifiable
        }
        bool IsRW() const
        {
            return (this->cFlags() & (NOTREAD | NOTWRITE)) == 0;    /// This funcion checks whether the vertex is both readable and modifiable
        }
        bool IsS() const
        {
            return (this->cFlags() & SELECTED) != 0;    ///  checks if the vertex is Selected
        }
        bool IsB() const
        {
            return (this->cFlags() & BORDER) != 0;    ///  checks if the vertex is a border one
        }
        bool IsV() const
        {
            return (this->cFlags() & VISITED) != 0;    ///  checks if the vertex Has been visited
        }


        /** Set the flag value
            @param flagp Valore da inserire nel flag
        */
        void SetFlags(int flagp)
        {
            this->Flags() = flagp;
        }

        /** Set the flag value
            @param flagp Valore da inserire nel flag
        */
        void ClearFlags()
        {
            this->Flags() = 0;
        }
        void SetD()
        {
            this->Flags() |= DELETED;    ///  deletes the vertex from the mesh
        }
        void ClearD()
        {
            this->Flags() &= (~DELETED);    ///  un-delete a vertex
        }
        void SetR()
        {
            this->Flags() &= (~NOTREAD);    ///  marks the vertex as readable
        }
        void ClearR()
        {
            this->Flags() |= NOTREAD;    ///  marks the vertex as not readable
        }
        void ClearW()
        {
            this->Flags() |= NOTWRITE;    ///  marks the vertex as writable
        }
        void SetW()
        {
            this->Flags() &= (~NOTWRITE);    ///  marks the vertex as not writable
        }
        void SetS()
        {
            this->Flags() |= SELECTED;    ///  select the vertex
        }
        void ClearS()
        {
            this->Flags() &= ~SELECTED;    /// Un-select a vertex
        }
        void SetB()
        {
            this->Flags() |= BORDER;
        }
        void ClearB()
        {
            this->Flags() &= ~BORDER;
        }
        void SetV()
        {
            this->Flags() |= VISITED;
        }
        void ClearV()
        {
            this->Flags() &= ~VISITED;
        }

        ///  Return the first bit that is not still used
        static int &FirstUnusedBitFlag()
        {
            static int b = USER0;
            return b;
        }

        /// Allocate a bit among the flags that can be used by user. It updates the FirstUnusedBitFlag.
        static inline int NewBitFlag()
        {
            int bitForTheUser = FirstUnusedBitFlag();
            FirstUnusedBitFlag() = FirstUnusedBitFlag() << 1;
            return bitForTheUser;
        }

        /// De-allocate a pre allocated bit. It updates the FirstUnusedBitFlag.
        // Note you must deallocate bit in the inverse order of the allocation (as in a stack)
        static inline bool DeleteBitFlag(int bitval)
        {
            if (FirstUnusedBitFlag() >> 1 == bitval)
            {
                FirstUnusedBitFlag() = FirstUnusedBitFlag() >> 1;
                return true;
            }
            assert(0);
            return false;
        }

        /// This function checks if the given user bit is true
        bool IsUserBit(int userBit)
        {
            return (this->Flags() & userBit) != 0;
        }

        /// This function set the given user bit
        void SetUserBit(int userBit)
        {
            this->Flags() |= userBit;
        }

        /// This function clear the given user bit
        void ClearUserBit(int userBit)
        {
            this->Flags() &= (~userBit);
        }

        template<class BoxType>
        void GetBBox(BoxType &bb) const
        {
            bb.Set(this->cP());
        }

    };


    template < class UserTypes,
        template <typename> class A = DefaultDeriver, template <typename> class B = DefaultDeriver,
        template <typename> class C = DefaultDeriver, template <typename> class D = DefaultDeriver,
        template <typename> class E = DefaultDeriver, template <typename> class F = DefaultDeriver >
    class Vertex : public VertexArityMax<UserTypes, A, B, C, D, E, F>
    {
    public:
        typedef AllTypes::AVertexType IAm;
        typedef UserTypes TypesPool;
    };


    /*------------------------------------------------------------------*/
    /*
    The base class of all the recusive definition chain. It is just a container of the typenames of the various simplexes.
    These typenames must be known form all the derived classes.
    */

    template <class UserTypes>
    class FaceTypeHolder : public UserTypes
    {
    public:

        template <class LeftF>
        void ImportData(const LeftF &) {}
        static void Name(std::vector<std::string> & /* name */) {}


        // prot
        inline int VN()  const
        {
            return 3;
        }
        inline int Prev(const int &i) const
        {
            return (i + (3 - 1)) % 3;
        }
        inline int Next(const int &i) const
        {
            return (i + 1) % 3;
        }
        inline void Alloc(const int &) {}
        inline void Dealloc() {}
    };

    namespace face
    {
        /** \addtogroup FaceComponentGroup
        @{
        */
        /*------------------------- EMPTY CORE COMPONENTS -----------------------------------------*/

        template <class T> class EmptyCore : public T
        {
        public:
            typedef tessellation::TexCoord2<float, 1> TexCoordType;

            inline typename T::CoordType cP(const int) const { assert(0);		static typename T::CoordType coord(0, 0, 0); return coord; }

            static bool HasMark() { return false; }

            TexCoordType &WT(const int)
            {
                static TexCoordType dummy_texture;
                assert(0);
                return dummy_texture;
            }
            TexCoordType const &cWT(const int) const
            {
                static TexCoordType dummy_texture;
                return dummy_texture;
            }

            typedef int FlagType;
            int &Flags()
            {
                static int dummyflags(0);
                assert(0);
                return dummyflags;
            }
            int cFlags() const
            {
                return 0;
            }
            static bool HasFlags()
            {
                return false;
            }

            typedef int VFAdjType;
            typename T::FacePointer &VFp(int)
            {
                static typename T::FacePointer fp = 0;
                assert(0);
                return fp;
            }
            typename T::FacePointer cVFp(int) const
            {
                static typename T::FacePointer fp = 0;
                assert(0);
                return fp;
            }
            typename T::FacePointer &FFp(int)
            {
                static typename T::FacePointer fp = 0;
                assert(0);
                return fp;
            }
            typename T::FacePointer cFFp(int) const
            {
                static typename T::FacePointer fp = 0;
                assert(0);
                return fp;
            }

            char &VFi(int)
            {
                static char z = 0;
                assert(0);
                return z;
            }
            char &FFi(int)
            {
                static char z = 0;
                assert(0);
                return z;
            }
            char cVFi(int) const
            {
                static char z = 0;
                assert(0);
                return z;
            }
            char cFFi(int) const
            {
                static char z = 0;
                assert(0);
                return z;
            }

            static bool HasVFAdjacency()
            {
                return false;
            }
            static bool HasFFAdjacency()
            {
                return false;
            }
            static bool HasFEAdjacency()
            {
                return false;
            }
            static bool HasFHAdjacency()
            {
                return false;
            }

            static bool HasPolyInfo() { return false; }

            bool IsVFInitialized(const int j) const { return  static_cast<const typename T::FaceType *>(this)->cVFi(j) != -1; }

            void VFClear(int j) {
                if (IsVFInitialized(j)) {
                    static_cast<typename T::FacePointer>(this)->VFp(j) = 0;
                    static_cast<typename T::FacePointer>(this)->VFi(j) = -1;
                }
            }
        };


        template <class VALUE_TYPE>
        class vector_ocf : public std::vector<VALUE_TYPE>
        {
            typedef std::vector<VALUE_TYPE> BaseType;
            typedef typename vector_ocf<VALUE_TYPE>::iterator ThisTypeIterator;

        public:
            vector_ocf() :std::vector<VALUE_TYPE>()
            {
                MarkEnabled = false;
                NormalEnabled = false;
                VFAdjacencyEnabled = false;
                FFAdjacencyEnabled = false;
            }

            // Auxiliary types to build internal vectors
            struct AdjTypePack {
                typename VALUE_TYPE::FacePointer _fp[3];
                char _zp[3];

                // Default constructor.
                // Needed because we need to know if adjacency is initialized or not
                // when resizing vectors and during an allocate face.
                AdjTypePack() {
                    _fp[0] = 0;
                    _fp[1] = 0;
                    _fp[2] = 0;
                }
            };


            ////////////////////////////////////////
            // All the standard methods of std::vector that can change the reallocation are
            // redefined in order to manage the additional data.
            void push_back(const VALUE_TYPE & v)
            {
                BaseType::push_back(v);
                BaseType::back()._ovp = this;
                if (MarkEnabled)        MV.push_back(0);
                if (NormalEnabled)      NV.push_back(typename VALUE_TYPE::NormalType());
                if (VFAdjacencyEnabled) AV.push_back(AdjTypePack());
                if (FFAdjacencyEnabled) AF.push_back(AdjTypePack());
            }
            void pop_back();

            virtual void resize(size_t _size)
            {
                size_t oldsize = BaseType::size();
                BaseType::resize(_size);
                if (oldsize < _size) {
                    ThisTypeIterator firstnew = BaseType::begin();
                    advance(firstnew, oldsize);
                    _updateOVP(firstnew, (*this).end());
                }
                if (MarkEnabled)        MV.resize(_size);
                if (VFAdjacencyEnabled) AV.resize(_size);
                if (FFAdjacencyEnabled) AF.resize(_size);
                if (NormalEnabled)      NV.resize(_size);
            }
            void reserve(size_t _size)
            {
                BaseType::reserve(_size);

                if (MarkEnabled)        MV.reserve(_size);
                if (NormalEnabled)      NV.reserve(_size);
                if (VFAdjacencyEnabled) AV.reserve(_size);
                if (FFAdjacencyEnabled) AF.reserve(_size);

                if (BaseType::empty()) return;

                ThisTypeIterator oldbegin = (*this).begin();
                if (oldbegin != (*this).begin()) _updateOVP((*this).begin(), (*this).end());
            }

            void _updateOVP(ThisTypeIterator lbegin, ThisTypeIterator lend)
            {
                ThisTypeIterator fi;
                //for(fi=(*this).begin();vi!=(*this).end();++vi)
                for (fi = lbegin; fi != lend; ++fi)
                    (*fi)._ovp = this;
            }

            bool IsMarkEnabled() const { return MarkEnabled; }
            void EnableMark() {
                assert(VALUE_TYPE::HasMarkOcf());
                MarkEnabled = true;
                MV.resize((*this).size(), 0);
            }

            void DisableMark() {
                assert(VALUE_TYPE::HasMarkOcf());
                MarkEnabled = false;
                MV.clear();
            }

            bool IsNormalEnabled() const { return NormalEnabled; }
            void EnableNormal() {
                assert(VALUE_TYPE::HasNormalOcf());
                NormalEnabled = true;
                NV.resize((*this).size());
            }

            void DisableNormal() {
                assert(VALUE_TYPE::HasNormalOcf());
                NormalEnabled = false;
                NV.clear();
            }

            bool IsVFAdjacencyEnabled() const { return VFAdjacencyEnabled; }
            void EnableVFAdjacency() {
                assert(VALUE_TYPE::HasVFAdjacencyOcf());
                VFAdjacencyEnabled = true;
                AV.resize((*this).size());
            }

            void DisableVFAdjacency() {
                assert(VALUE_TYPE::HasVFAdjacencyOcf());
                VFAdjacencyEnabled = false;
                AV.clear();
            }


            bool IsFFAdjacencyEnabled() const { return FFAdjacencyEnabled; }
            void EnableFFAdjacency() {
                assert(VALUE_TYPE::HasFFAdjacencyOcf());
                FFAdjacencyEnabled = true;
                AF.resize((*this).size());
            }

            void DisableFFAdjacency() {
                assert(VALUE_TYPE::HasFFAdjacencyOcf());
                FFAdjacencyEnabled = false;
                AF.clear();
            }

        public:
            std::vector<typename VALUE_TYPE::NormalType> NV;
            std::vector<int> MV;
            std::vector<struct AdjTypePack> AV;
            std::vector<struct AdjTypePack> AF;

            bool MarkEnabled;
            bool NormalEnabled;
            bool VFAdjacencyEnabled;
            bool FFAdjacencyEnabled;
        }; // end class vector_ocf


        ///*-------------------------- InfoOpt  ----------------------------------*/
        template < class T> class InfoOcf : public T {
        public:
            // You should never ever try to copy a vertex that has OCF stuff.
            // use ImportData function.
            inline InfoOcf &operator=(const InfoOcf & /*other*/) {
                assert(0); return *this;
            }


            vector_ocf<typename T::FaceType> &Base() const
            {
                return *_ovp;
            }

            template <class RightFaceType>
            void ImportData(const RightFaceType & rightF) { T::ImportData(rightF); }

            static bool HasMarkOcf() { return false; }
            static bool HasNormalOcf() { return false; }
            static bool HasFFAdjacencyOcf() { return false; }
            static bool HasVFAdjacencyOcf() { return false; }


            inline size_t Index() const {
                typename T::FaceType const *tp = static_cast<typename T::FaceType const *>(this);
                size_t tt2 = tp - &*(_ovp->begin());
                return tt2;
            }
        public:
            // ovp Optional Vector Pointer
            // Pointer to the base vector where each face element is stored.
            // used to access to the vectors of the other optional members.
            vector_ocf<typename T::FaceType> *_ovp;
        };

        template < class FaceType >
        bool FaceVectorHasVFAdjacency(const face::vector_ocf<FaceType> &fv)
        {
            if (FaceType::HasVFAdjacencyOcf()) return fv.IsVFAdjacencyEnabled();
            else return FaceType::HasVFAdjacency();
        }
        template < class FaceType >
        bool FaceVectorHasFFAdjacency(const face::vector_ocf<FaceType> &fv)
        {
            if (FaceType::HasFFAdjacencyOcf()) return fv.IsFFAdjacencyEnabled();
            else return FaceType::HasFFAdjacency();
        }

        /*----------------------------- VFADJ ------------------------------*/
        template <class T> class VFAdjOcf : public T {
        public:
            typename T::FacePointer &VFp(const int j) {
                assert((*this).Base().VFAdjacencyEnabled);
                return (*this).Base().AV[(*this).Index()]._fp[j];
            }

            typename T::FacePointer cVFp(const int j) const {
                if (!(*this).Base().VFAdjacencyEnabled) return 0;
                else return (*this).Base().AV[(*this).Index()]._fp[j];
            }

            char &VFi(const int j) {
                assert((*this).Base().VFAdjacencyEnabled);
                return (*this).Base().AV[(*this).Index()]._zp[j];
            }

            char cVFi(const int j) const {
                assert((*this).Base().VFAdjacencyEnabled);
                return (*this).Base().AV[(*this).Index()]._zp[j];
            }

            template <class RightFaceType>
            void ImportData(const RightFaceType & rightF) {
                T::ImportData(rightF);
            }
            static bool HasVFAdjacency() { return true; }
            static bool HasVFAdjacencyOcf() { return true; }

        private:
        };

        /*----------------------------- FFADJ ------------------------------*/
        template <class T> class FFAdjOcf : public T {
        public:
            typename T::FacePointer &FFp(const int j) {
                assert((*this).Base().FFAdjacencyEnabled);
                return (*this).Base().AF[(*this).Index()]._fp[j];
            }

            typename T::FacePointer cFFp(const int j) const {
                if (!(*this).Base().FFAdjacencyEnabled) return 0;
                else return (*this).Base().AF[(*this).Index()]._fp[j];
            }

            char &FFi(const int j) {
                assert((*this).Base().FFAdjacencyEnabled);
                return (*this).Base().AF[(*this).Index()]._zp[j];
            }
            char cFFi(const int j) const {
                assert((*this).Base().FFAdjacencyEnabled);
                return (*this).Base().AF[(*this).Index()]._zp[j];
            }

            typename T::FacePointer  &FFp1(const int j) { return FFp((j + 1) % 3); }
            typename T::FacePointer  &FFp2(const int j) { return FFp((j + 2) % 3); }
            typename T::FacePointer  cFFp1(const int j) const { return FFp((j + 1) % 3); }
            typename T::FacePointer  cFFp2(const int j) const { return FFp((j + 2) % 3); }

            typename T::FacePointer  &Neigh(const int j) { return FFp(j); }
            typename T::FacePointer  cNeigh(const int j) const { return cFFp(j); }
            unsigned int SizeNeigh() { return 3; }

            template <class RightFaceType>
            void ImportData(const RightFaceType & rightF) {
                T::ImportData(rightF);
            }
            static bool HasFFAdjacency() { return true; }
            static bool HasFFAdjacencyOcf() { return true; }
        };

        template <class A, class T> class NormalAbs : public T {
        public:
            typedef A NormalType;
            inline NormalType &N() { return _norm; }
            inline NormalType cN() const { return _norm; }

            inline void Alloc(const int & ns) { T::Alloc(ns); }
            inline void Dealloc() { T::Dealloc(); }
            static bool HasNormal() { return true; }
            static void Name(std::vector<std::string> & name) { name.push_back(std::string("NormalAbs")); T::Name(name); }

        private:
            NormalType _norm;
        };

        template <class T> class Normal3f : public NormalAbs<tessellation::Point3f, T> {
        public:  static void Name(std::vector<std::string> & name) { name.push_back(std::string("Normal3f")); T::Name(name); }
        };
        

        /*------------------------- BitFlags -----------------------------------------*/
        /*! \brief \em Component: Per face \b Flags

        This component stores a 32 bit array of bit flags. These bit flags are used for keeping track of selection, deletion, visiting etc. \sa \ref flags for more details on common uses of flags.
        */
        template <class T> class BitFlags : public T {
        public:
            BitFlags() :_flags(0) {}
            typedef int FlagType;
            int &Flags() { return _flags; }
            int cFlags() const { return _flags; }
            template <class RightValueType>
            void ImportData(const RightValueType & rightF) {
                if (RightValueType::HasFlags())
                    Flags() = rightF.cFlags();
                T::ImportData(rightF);
            }
            inline void Alloc(const int & ns) { T::Alloc(ns); }
            inline void Dealloc() { T::Dealloc(); }
            static bool HasFlags() { return true; }
            static void Name(std::vector<std::string> & name) { name.push_back(std::string("BitFlags")); T::Name(name); }

        private:
            int  _flags;
        };

        /*-------------------------- VertexRef ----------------------------------------*/
        /*! \brief The references to the vertexes of a triangular face
        *
        * Stored as three pointers to the VertexType
        */


        template <class T> class VertexRef : public T
        {
        public:
            typedef typename T::VertexType::CoordType CoordType;
            typedef typename T::VertexType::ScalarType ScalarType;

            VertexRef()
            {
                v[0] = 0;
                v[1] = 0;
                v[2] = 0;
            }

            inline CoordType cP(const int j) const { assert(j >= 0 && j < 3);		return v[j]->cP(); }

            inline typename T::VertexType*&V(const int j)
            {
                assert(j >= 0 && j < 3);    /// \brief The pointer to the i-th vertex
                return v[j];
            }
            inline typename T::VertexType *cV(const int j) const
            {
                assert(j >= 0 && j < 3);
                return v[j];
            }

            inline typename T::VertexType *&V0(const int j)
            {
                return V(j);    /** \brief Return the pointer to the j-th vertex of the face. */
            }
            inline typename T::VertexType *&V1(const int j)
            {
                return V((j + 1) % 3);    /** \brief Return the pointer to the ((j+1)%3)-th vertex of the face. */
            }
            inline typename T::VertexType *&V2(const int j)
            {
                return V((j + 2) % 3);    /** \brief Return the pointer to the ((j+2)%3)-th vertex of the face. */
            }


        private:
            typename T::VertexType *v[3];
        };

        template <class FaceType>
        class Pos
        {
        public:

            /// The vertex type
            typedef typename FaceType::VertexType VertexType;
            ///The Pos type
            typedef Pos<FaceType> PosType;
            /// The scalar type
            typedef typename VertexType::ScalarType ScalarType;

            /// Pointer to the face of the half-edge
            typename FaceType::FaceType *f;
            /// Index of the edge
            int z;
            /// Pointer to the vertex
            VertexType *v;

            /// Default constructor
            Pos() : f(0), z(-1), v(0) {}
            /// Constructor which associates the half-edge element with a face, its edge and its vertex
            /// \note that the input must be consistent, e.g. it should hold that \c vp==fp->V0(zp) or \c vp==fp->V1(zp)
            Pos(FaceType *const fp, int const zp, VertexType *const vp)
            {
                f = fp;
                z = zp;
                v = vp;
                assert((vp == fp->V0(zp)) || (vp == fp->V1(zp)));
            }
            Pos(FaceType *const fp, int const zp)
            {
                f = fp;
                z = zp;
                v = f->V(zp);
            }
            Pos(FaceType *const fp, VertexType *const vp)
            {
                f = fp;
                v = vp;
                for (int i = 0; i < f->VN(); ++i)
                    if (f->V(i) == v)
                    {
                        z = f->Prev(i);
                        break;
                    }
            }

            // Official Access functions functions
            VertexType  *&V()
            {
                return v;
            }
            int          &E()
            {
                return z;
            }
            FaceType    *&F()
            {
                return f;
            }

            VertexType *V() const
            {
                return v;
            }
            int          E() const
            {
                return z;
            }
            FaceType    *F() const
            {
                return f;
            }

            /// Operator to compare two half-edge
            inline bool operator == (PosType const &p) const
            {
                return (f == p.f && z == p.z && v == p.v);
            }

            /// Operator to compare two half-edge
            inline bool operator != (PosType const &p) const
            {
                return (f != p.f || z != p.z || v != p.v);
            }


            //Cambia Faccia lungo z
            // e' uguale a FlipF solo che funziona anche per non manifold.
            /// Change face via z
            void NextF()
            {
                FaceType *t = f;
                f = t->FFp(z);
                z = t->FFi(z);
            }

            // Paolo Cignoni 19/6/99
            // Si muove sulla faccia adiacente a f, lungo uno spigolo che
            // NON e' j, e che e' adiacente a v
            // in questo modo si scandiscono tutte le facce incidenti in un
            // vertice f facendo Next() finche' non si ritorna all'inizio
            // Nota che sul bordo rimbalza, cioe' se lo spigolo !=j e' di bordo
            // restituisce sempre la faccia f ma con nj che e' il nuovo spigolo di bordo
            // vecchi parametri:     	FaceType * & f, VertexType * v, int & j

            /// It moves on the adjacent face incident to v, via a different edge that j
            void NextE()
            {
                assert(f->V(z) == v || f->V(f->Next(z)) == v); // L'edge j deve contenere v
                FlipE();
                FlipF();
                assert(f->V(z) == v || f->V(f->Next(z)) == v);
            }
            // Cambia edge mantenendo la stessa faccia e lo stesso vertice
            /// Changes edge maintaining the same face and the same vertex
            void FlipE()
            {
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V((z + 0) % f->VN()) == v));
                if (f->V(f->Next(z)) == v) z = f->Next(z);
                else z = f->Prev(z);
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V((z)) == v));
            }

            // Cambia Faccia mantenendo lo stesso vertice e lo stesso edge
            // Vale che he.flipf.flipf= he
            // Se l'he e' di bordo he.flipf()==he
            // Si puo' usare SOLO se l'edge e' 2manifold altrimenti
            // si deve usare nextf

            /// Changes face maintaining the same vertex and the same edge
            void FlipF()
            {
                assert(f->FFp(z)->FFp(f->FFi(z)) == f);  // two manifoldness check
                // Check that pos vertex is one of the current z-th edge and it is different from the vert opposite to the edge.
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V((z)) == v));
                FaceType *nf = f->FFp(z);
                int nz = f->FFi(z);
                assert(nf->V(nf->Prev(nz)) != v && (nf->V(nf->Next(nz)) == v || nf->V((nz)) == v));
                f = nf;
                z = nz;
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));
            }

            /// Changes vertex maintaining the same face and the same edge
            void FlipV()
            {
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));

                if (f->V(f->Next(z)) == v)
                    v = f->V(z);
                else
                    v = f->V(f->Next(z));

                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));
            }

            /// return the vertex that it should have if we make FlipV;
            VertexType *VFlip() const
            {
                assert(f->cV(f->Prev(z)) != v && (f->cV(f->Next(z)) == v || f->cV(z) == v));
                if (f->cV(f->Next(z)) == v)	return f->cV(z);
                else			return f->cV(f->Next(z));
            }

            /// return the face that it should have if we make FlipF;
            FaceType *FFlip() const
            {
                //        assert( f->FFp(z)->FFp(f->FFi(z))==f );
                //        assert(f->V(f->Prev(z))!=v);
                //        assert(f->V(f->Next(z))==v || f->V((z+0)%f->VN())==v);
                FaceType *nf = f->FFp(z);
                return nf;
            }


            // Trova il prossimo half-edge di bordo (nhe)
            // tale che
            // --nhe.f adiacente per vertice a he.f
            // --nhe.v adiacente per edge di bordo a he.v
            // l'idea e' che se he e' un half edge di bordo
            // si puo scorrere tutto un bordo facendo
            //
            //		hei=he;
            //		do
            //			hei.Nextb()
            //		while(hei!=he);

            /// Finds the next half-edge border
            void NextB()
            {
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));
                assert(f->FFp(z) == f); // f is border along j
                // Si deve cambiare faccia intorno allo stesso vertice v
                //finche' non si trova una faccia di bordo.
                do
                    NextE();
                while (!IsBorder());

                // L'edge j e' di bordo e deve contenere v
                assert(IsBorder() && (f->V(z) == v || f->V(f->Next(z)) == v));

                FlipV();
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));
                assert(f->FFp(z) == f); // f is border along j
            }

            /// Checks if the half-edge is of border
            bool IsBorder()const
            {
                return face::IsBorder(*f, z);
            }

        };

        /** Return a boolean that indicate if the j-th edge of the face is a border.
        @param j Index of the edge
        @return true if j is an edge of border, false otherwise
        */
        template <class FaceType>
        inline bool IsBorder(FaceType const &f, const int j)
        {
            if (FaceType::HasFFAdjacency())
                return f.cFFp(j) == &f;
            //return f.IsBorder(j);

            assert(0);
            return true;
        }

        ///*-------------------------- MARK  ----------------------------------*/
        template <class T> class MarkOcf : public T {
        public:
            inline int &IMark() {
                assert((*this).Base().MarkEnabled);
                return (*this).Base().MV[(*this).Index()];
            }
            inline int cIMark() const {
                assert((*this).Base().MarkEnabled);
                return (*this).Base().MV[(*this).Index()];
            }

            template <class RightFaceType>
            void ImportData(const RightFaceType & rightF) {
                if ((*this).IsMarkEnabled() && rightF.IsMarkEnabled())
                    IMark() = rightF.cIMark();
                T::ImportData(rightF);
            }
            inline bool IsMarkEnabled()          const { return this->Base().IsMarkEnabled(); }
            static bool HasMark() { return true; }
            static bool HasMarkOcf() { return true; }
            inline void InitIMark() { IMark() = 0; }
        };

        /** @} */   // End Doxygen FaceComponentGroup

        template <class S>
        class PointDistanceBaseFunctor {
        public:
            typedef S ScalarType;
            typedef Point3<ScalarType> QueryType;

            static inline const Point3<ScalarType> & Pos(const Point3<ScalarType> & qt) { return qt; }
            template <class FACETYPE, class SCALARTYPE>
            inline bool operator () (const FACETYPE & f, const Point3<SCALARTYPE> & p, SCALARTYPE & minDist, Point3<SCALARTYPE> & q) {
                const Point3<typename FACETYPE::ScalarType> fp = Point3<typename FACETYPE::ScalarType>::Construct(p);
                Point3<typename FACETYPE::ScalarType> fq;
                typename FACETYPE::ScalarType md = (typename FACETYPE::ScalarType)(minDist);
                const bool ret = PointDistanceBase(f, fp, md, fq);
                minDist = (SCALARTYPE)(md);
                q = Point3<SCALARTYPE>::Construct(fq);
                return (ret);
            }
        };

        /// BASIC VERSION of the Point-face distance that does not require the EdgePlane Additional data.
        /// Given a face and a point, returns the closest point of the face to p.

        template <class FaceType>
        bool PointDistanceBase(
            const FaceType &f,																		/// the face to be tested
            const tessellation::Point3<typename FaceType::ScalarType> & q, /// the point tested
            typename FaceType::ScalarType & dist,                 /// bailout distance. It must be initialized with the max admittable value.
            tessellation::Point3<typename FaceType::ScalarType> & p)
        {
            typedef typename FaceType::ScalarType ScalarType;

            if (f.cN() == Point3<ScalarType>(0, 0, 0)) // to correctly manage the case of degenerate triangles we consider them as segments.
            {
                Box3<ScalarType> bb;
                f.GetBBox(bb);
                Segment3<ScalarType> degenTri(bb.min, bb.max);
                Point3<ScalarType> closest;
                ScalarType d;
                if (bb.Diag() > 0)
                    tessellation::SegmentPointDistance<ScalarType>(degenTri, q, closest, d);
                else // very degenerate triangle (just a point)
                {
                    closest = bb.min;
                    d = Distance(q, closest);
                }
                if (d > dist) return false;
                dist = d;
                p = closest;
                assert(!math::IsNAN(dist));
                return true;
            }

            Plane3<ScalarType, true> fPlane;
            fPlane.Init(f.cP(0), f.cN());
            const ScalarType EPS = ScalarType(0.000001);
            ScalarType b, b0, b1, b2;
            // Calcolo distanza punto piano
            ScalarType d = SignedDistancePlanePoint(fPlane, q);
            if (d > dist || d < -dist)			// Risultato peggiore: niente di fatto
                return false;

            // Projection of query point onto the triangle plane
            p = q - fPlane.Direction()*d;

            Point3<ScalarType> fEdge[3];
            fEdge[0] = f.cP(1); fEdge[0] -= f.cP(0);
            fEdge[1] = f.cP(2); fEdge[1] -= f.cP(1);
            fEdge[2] = f.cP(0); fEdge[2] -= f.cP(2);


            int bestAxis;
            if (fabs(f.cN()[0]) > fabs(f.cN()[1]))
            {
                if (fabs(f.cN()[0]) > fabs(f.cN()[2])) bestAxis = 0;
                else bestAxis = 2;
            }
            else {
                if (fabs(f.cN()[1]) > fabs(f.cN()[2])) bestAxis = 1; /* 1 > 0 ? 2 */
                else bestAxis = 2; /* 2 > 1 ? 2 */
            }

            ScalarType scaleFactor;

            switch (bestAxis)
            {
            case 0:  /************* X AXIS **************/
                scaleFactor = 1 / fPlane.Direction()[0];
                fEdge[0] *= scaleFactor; fEdge[1] *= scaleFactor; fEdge[2] *= scaleFactor;

                b0 = fEdge[1][1] * (p[2] - f.cP(1)[2]) - fEdge[1][2] * (p[1] - f.cP(1)[1]);
                if (b0 <= 0)
                {
                    b0 = PSDist(q, f.cV(1)->cP(), f.cV(2)->cP(), p);
                    if (dist > b0) { dist = b0; return true; }
                    else return false;
                }
                b1 = fEdge[2][1] * (p[2] - f.cP(2)[2]) - fEdge[2][2] * (p[1] - f.cP(2)[1]);
                if (b1 <= 0)
                {
                    b1 = PSDist(q, f.cV(2)->cP(), f.cV(0)->cP(), p);
                    if (dist > b1) { dist = b1; return true; }
                    else return false;
                }
                b2 = fEdge[0][1] * (p[2] - f.cP(0)[2]) - fEdge[0][2] * (p[1] - f.cP(0)[1]);
                if (b2 <= 0)
                {
                    b2 = PSDist(q, f.cV(0)->cP(), f.cV(1)->cP(), p);
                    if (dist > b2) { dist = b2; return true; }
                    else return false;
                }

                if ((b = math::Min<ScalarType>(b0, b1, b2)) < EPS*DoubleArea(f))
                {
                    ScalarType bt;
                    if (b == b0) 	    bt = PSDist(q, f.cV(1)->cP(), f.cV(2)->cP(), p);
                    else if (b == b1) 	bt = PSDist(q, f.cV(2)->cP(), f.cV(0)->cP(), p);
                    else { assert(b == b2); bt = PSDist(q, f.cV(0)->cP(), f.cV(1)->cP(), p); }
                    //printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
                    if (dist > bt) { dist = bt; return true; }
                    else return false;
                }
                break;

            case 1:  /************* Y AXIS **************/
                scaleFactor = 1 / fPlane.Direction()[1];
                fEdge[0] *= scaleFactor; fEdge[1] *= scaleFactor; fEdge[2] *= scaleFactor;

                b0 = fEdge[1][2] * (p[0] - f.cP(1)[0]) - fEdge[1][0] * (p[2] - f.cP(1)[2]);
                if (b0 <= 0)
                {
                    b0 = PSDist(q, f.cV(1)->cP(), f.cV(2)->cP(), p);
                    if (dist > b0) { dist = b0; return true; }
                    else return false;
                }
                b1 = fEdge[2][2] * (p[0] - f.cP(2)[0]) - fEdge[2][0] * (p[2] - f.cP(2)[2]);
                if (b1 <= 0)
                {
                    b1 = PSDist(q, f.cV(2)->cP(), f.cV(0)->cP(), p);
                    if (dist > b1) { dist = b1; return true; }
                    else return false;
                }
                b2 = fEdge[0][2] * (p[0] - f.cP(0)[0]) - fEdge[0][0] * (p[2] - f.cP(0)[2]);
                if (b2 <= 0)
                {
                    b2 = PSDist(q, f.cV(0)->cP(), f.cV(1)->cP(), p);
                    if (dist > b2) { dist = b2; return true; }
                    else return false;
                }
                if ((b = math::Min<ScalarType>(b0, b1, b2)) < EPS*DoubleArea(f))
                {
                    ScalarType bt;
                    if (b == b0) 	    bt = PSDist(q, f.cV(1)->cP(), f.cV(2)->cP(), p);
                    else if (b == b1) 	bt = PSDist(q, f.cV(2)->cP(), f.cV(0)->cP(), p);
                    else { assert(b == b2); bt = PSDist(q, f.cV(0)->cP(), f.cV(1)->cP(), p); }
                    //printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);
                    if (dist > bt) { dist = bt; return true; }
                    else return false;
                }
                break;

            case 2:  /************* Z AXIS **************/
                scaleFactor = 1 / fPlane.Direction()[2];
                fEdge[0] *= scaleFactor; fEdge[1] *= scaleFactor; fEdge[2] *= scaleFactor;

                b0 = fEdge[1][0] * (p[1] - f.cP(1)[1]) - fEdge[1][1] * (p[0] - f.cP(1)[0]);
                if (b0 <= 0)
                {
                    b0 = PSDist(q, f.cV(1)->cP(), f.cV(2)->cP(), p);
                    if (dist > b0) { dist = b0; return true; }
                    else return false;
                }
                b1 = fEdge[2][0] * (p[1] - f.cP(2)[1]) - fEdge[2][1] * (p[0] - f.cP(2)[0]);
                if (b1 <= 0)
                {
                    b1 = PSDist(q, f.cV(2)->cP(), f.cV(0)->cP(), p);
                    if (dist > b1) { dist = b1; return true; }
                    else return false;
                }
                b2 = fEdge[0][0] * (p[1] - f.cP(0)[1]) - fEdge[0][1] * (p[0] - f.cP(0)[0]);
                if (b2 <= 0)
                {
                    b2 = PSDist(q, f.cV(0)->cP(), f.cV(1)->cP(), p);
                    if (dist > b2) { dist = b2; return true; }
                    else return false;
                }
                if ((b = math::Min<ScalarType>(b0, b1, b2)) < EPS*DoubleArea(f))
                {
                    ScalarType bt;
                    if (b == b0) 	    bt = PSDist(q, f.cV(1)->cP(), f.cV(2)->cP(), p);
                    else if (b == b1) 	bt = PSDist(q, f.cV(2)->cP(), f.cV(0)->cP(), p);
                    else { assert(b == b2); bt = PSDist(q, f.cV(0)->cP(), f.cV(1)->cP(), p); }
                    //printf("Warning area:%g %g %g %g thr:%g bt:%g\n",Area(), b0,b1,b2,EPSILON*Area(),bt);

                    if (dist > bt) { dist = bt; return true; }
                    else return false;
                }
                break;
            default: assert(0); // if you get this assert it means that you forgot to set the required UpdateFlags<MeshType>::FaceProjection(m);

            }

            dist = ScalarType(fabs(d));
            //dist = Distance(p,q);
            return true;
        }

    } // end namespace face

    /* The base class form which we start to add our components.
    it has the empty definition for all the standard members (coords, color flags)
    Note:
    in order to avoid both virtual classes and ambiguous definitions all
    the subsequent overrides must be done in a sequence of derivation.

    In other words we cannot derive and add in a single derivation step
    (with multiple ancestor), both the real (non-empty) normal and color but
    we have to build the type a step a time (deriving from a single ancestor at a time).


    */
    template <class UserTypes>
    class FaceBase : public
        face::EmptyCore< FaceTypeHolder <UserTypes> >
    {
    };


    /* The Real Big Face class;

    The class __FaceArityMax__ is the one that is the Last to be derived,
    and therefore is the only one to know the real members
    (after the many overrides) so all the functions with common behaviour
    using the members defined in the various Empty/nonEmpty component classes
    MUST be defined here.

    I.e. IsD() that uses the overridden Flags() member must be defined here.

    */

    template < class UserTypes,
        template <typename> class A, template <typename> class B, template <typename> class C,
        template <typename> class D, template <typename> class E,
        template <typename> class F, template <typename> class G>
    class FaceArityMax : public  Arity7<FaceBase<UserTypes>, A, B, C, D, E, F, G>
    {

    public:
        typedef typename  FaceArityMax::ScalarType ScalarType;
        // ----- Flags stuff -----

        enum
        {
            DELETED = 0x00000001,		// Face is deleted from the mesh
            NOTREAD = 0x00000002,		// Face of the mesh is not readable
            NOTWRITE = 0x00000004,		// Face of the mesh is not writable
            VISITED = 0x00000010,		// Face has been visited. Usualy this is a per-algorithm used bit.
            SELECTED = 0x00000020,		// Face is selected. Algorithms should try to work only on selected face (if explicitly requested)
            // Border _flags, it is assumed that BORDERi = BORDER0<<i
            BORDER0 = 0x00000040,
            BORDER1 = 0x00000080,
            BORDER2 = 0x00000100,
            BORDER012 = BORDER0 | BORDER1 | BORDER2,
 

            // Faux edges. (semantics: when a mesh is polygonal, edges which are inside a polygonal face are "faux"
            FAUX0 = 0x00040000,
            FAUX1 = 0x00080000,
            FAUX2 = 0x00100000,
            FAUX012 = FAUX0 | FAUX1 | FAUX2,
            // First user bit
            USER0 = 0x00200000
        };


        ///  checks if the Face is deleted
        bool IsD() const
        {
            return (this->cFlags() & DELETED) != 0;
        }
        ///  checks if the Face is readable
        bool IsR() const
        {
            return (this->cFlags() & NOTREAD) == 0;
        }
        ///  checks if the Face is modifiable
        bool IsW() const
        {
            return (this->cFlags() & NOTWRITE) == 0;
        }
        /// This funcion checks whether the Face is both readable and modifiable
        bool IsRW() const
        {
            return (this->cFlags() & (NOTREAD | NOTWRITE)) == 0;
        }
        ///  checks if the Face is Modified
        bool IsS() const
        {
            return (this->cFlags() & SELECTED) != 0;
        }
        ///  checks if the Face is Modified
        bool IsV() const
        {
            return (this->cFlags() & VISITED) != 0;
        }

        /** Set the flag value
        @param flagp Valore da inserire nel flag
        */
        void SetFlags(int flagp)
        {
            this->Flags() = flagp;
        }

        /** Set the flag value
        @param flagp Valore da inserire nel flag
        */
        void ClearFlags()
        {
            this->Flags() = 0;
        }
        
        ///  select the Face
        void SetS()
        {
            this->Flags() |= SELECTED;
        }
        /// Un-select a Face
        void ClearS()
        {
            this->Flags() &= ~SELECTED;
        }
        ///  select the Face
        void SetV()
        {
            this->Flags() |= VISITED;
        }
        /// Un-select a Face
        void ClearV()
        {
            this->Flags() &= ~VISITED;
        }

        /// This function checks if the face is selected
        bool IsB(int i) const
        {
            return (this->cFlags() & (BORDER0 << i)) != 0;
        }
        /// This function select the face
        void SetB(int i)
        {
            this->Flags() |= (BORDER0 << i);
        }
        /// This funcion execute the inverse operation of SetS()
        void ClearB(int i)
        {
            this->Flags() &= (~(BORDER0 << i));
        }

        /// This function checks if a given side of the face is a feature/internal edge
        /// it is used by some importer to mark internal
        /// edges of polygonal faces that have been triangulated
        bool IsF(int i) const
        {
            return (this->cFlags() & (FAUX0 << i)) != 0;
        }
        bool IsAnyF() const
        {
            return (this->cFlags() & (FAUX0 | FAUX1 | FAUX2)) != 0;
        }
        /// This function select the face
        void SetF(int i)
        {
            this->Flags() |= (FAUX0 << i);
        }
        /// This funcion execute the inverse operation of SetS()
        void ClearF(int i)
        {
            this->Flags() &= (~(FAUX0 << i));
        }
        void ClearAllF()
        {
            this->Flags() &= (~(FAUX0 | FAUX1 | FAUX2));
        }

        void GetBBox(Box3<ScalarType> &bb) const
        {
            if (this->IsD())
            {
                bb.SetNull();
                return;
            }
            bb.Set(this->cP(0));
            bb.Add(this->cP(1));
            bb.Add(this->cP(2));
        }
    };



    template < class UserTypes,
        template <typename> class A = DefaultDeriver, template <typename> class B = DefaultDeriver,
        template <typename> class C = DefaultDeriver, template <typename> class D = DefaultDeriver,
        template <typename> class E = DefaultDeriver, template <typename> class F = DefaultDeriver,
        template <typename> class G = DefaultDeriver>
    class Face : public FaceArityMax<UserTypes, A, B, C, D, E, F, G>
    {
    public:
        typedef AllTypes::AFaceType IAm;
        typedef UserTypes TypesPool;
    };

}
