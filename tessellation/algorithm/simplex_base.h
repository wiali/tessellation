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

        inline short      &n()
        {
            return _n[0];
        }
        inline short n() const
        {
            return _n[0];
        }

        inline Point2<T> &t()
        {
            return _t[0];
        }        
        inline Point2<T> t() const
        {
            return _t[0];
        }

        inline const PointType &P() const { return _t[0]; }
        inline PointType &P() { return _t[0]; }

        inline const PointType &P(const int i) const { assert(i > 0 && i < NMAX); return _t[i]; }
        inline PointType &P(const int i) { assert(i > 0 && i < NMAX); return _t[i]; }

        inline T & u() { return _t[0][0]; }
        inline T & v() { return _t[0][1]; }
        inline const T & u() const { return _t[0][0]; }
        inline const T & v() const { return _t[0][1]; }

        enum { n_coords = NMAX };
    };

    typedef TexCoord2<float>  TexCoord2f;

    namespace vertex
    {
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
            void ImportData(const RightValueType& /*rVert*/) {}

            static bool HasVEAdjacency() { return false; }
        };

        template <class A, class T> class Coord : public T
        {
        public:
            typedef A CoordType;
            typedef typename A::ScalarType      ScalarType;

            inline const CoordType &P() const
            {
                return _coord;
            }

            inline       CoordType &P()
            {
                return _coord;
            }

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

        template <class T> class Coord3f : public Coord<Point3f, T> {
        public:	static void Name(std::vector<std::string> & name) { name.push_back(std::string("Coord3f")); T::Name(name); }
        };

        template <class A, class T> class Normal : public T
        {
        public:
            typedef A NormalType;

            inline const NormalType &N() const
            {
                return _norm;
            }

            inline       NormalType &N()
            {
                return _norm;
            }

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

        template <class T> class Normal3f : public Normal<Point3f, T> {
        public:	static void Name(std::vector<std::string> & name) { name.push_back(std::string("Normal3f")); T::Name(name); }
        };

        template <class A, class TT> class TexCoord : public TT
        {
        public:
            typedef A TexCoordType;

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

        template <class TT> class TexCoord2f : public TexCoord<TexCoord2<float, 1>, TT>
        {
        public:
            static void Name(std::vector<std::string> &name)
            {
                name.push_back(std::string("TexCoord2f"));
                TT::Name(name);
            }
        };

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

        template <class A, class TT> class TexCoordOcf : public TT {
        public:
            typedef A TexCoordType;
            const TexCoordType &T() const { assert((*this).Base().TexCoordEnabled); return (*this).Base().TV[(*this).Index()]; }
            TexCoordType &T() { assert((*this).Base().TexCoordEnabled); return (*this).Base().TV[(*this).Index()]; }
            TexCoordType cT() const { assert((*this).Base().TexCoordEnabled); return (*this).Base().TV[(*this).Index()]; }
            template < class RightVertexType>
            void ImportData(const RightVertexType & rightV)
            {
                if ((*this).IsTexCoordEnabled() && rightV.IsTexCoordEnabled())
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

        template < class T> class InfoOcf : public T {
        public:
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
        };
    }

    template < class UserTypes,
        template <typename> class A, template <typename> class B,
        template <typename> class C, template <typename> class D,
        template <typename> class E, template <typename> class F>
    class VertexArityMax : public Arity6<vertex::EmptyCore<UserTypes>, A, B, C, D, E, F>
    {
    public:
        enum
        {
            DELETED = 0x0001,	
            NOTREAD = 0x0002,	
            NOTWRITE = 0x0004,	
            MODIFIED = 0x0008,	
            VISITED = 0x0010,	
            SELECTED = 0x0020,	
            BORDER = 0x0100,    
            USER0 = 0x0200		
        };

        bool IsD() const
        {
            return (this->cFlags() & DELETED) != 0;
        }
        bool IsR() const
        {
            return (this->cFlags() & NOTREAD) == 0;
        }
        bool IsW() const
        {
            return (this->cFlags() & NOTWRITE) == 0;
        }
        bool IsRW() const
        {
            return (this->cFlags() & (NOTREAD | NOTWRITE)) == 0;
        }
        bool IsS() const
        {
            return (this->cFlags() & SELECTED) != 0;
        }
        bool IsB() const
        {
            return (this->cFlags() & BORDER) != 0;
        }
        bool IsV() const
        {
            return (this->cFlags() & VISITED) != 0;
        }

        void SetFlags(int flagp)
        {
            this->Flags() = flagp;
        }

        void ClearFlags()
        {
            this->Flags() = 0;
        }
        void SetD()
        {
            this->Flags() |= DELETED;
        }
        void ClearD()
        {
            this->Flags() &= (~DELETED);
        }
        void SetR()
        {
            this->Flags() &= (~NOTREAD);
        }
        void ClearR()
        {
            this->Flags() |= NOTREAD;
        }
        void ClearW()
        {
            this->Flags() |= NOTWRITE;
        }
        void SetW()
        {
            this->Flags() &= (~NOTWRITE);
        }
        void SetS()
        {
            this->Flags() |= SELECTED;
        }
        void ClearS()
        {
            this->Flags() &= ~SELECTED;
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

        static int &FirstUnusedBitFlag()
        {
            static int b = USER0;
            return b;
        }

        static inline int NewBitFlag()
        {
            int bitForTheUser = FirstUnusedBitFlag();
            FirstUnusedBitFlag() = FirstUnusedBitFlag() << 1;
            return bitForTheUser;
        }

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

        bool IsUserBit(int userBit)
        {
            return (this->Flags() & userBit) != 0;
        }

        void SetUserBit(int userBit)
        {
            this->Flags() |= userBit;
        }

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

    template <class UserTypes>
    class FaceTypeHolder : public UserTypes
    {
    public:

        template <class LeftF>
        void ImportData(const LeftF &) {}
        static void Name(std::vector<std::string> & /* name */) {}

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

            struct AdjTypePack {
                typename VALUE_TYPE::FacePointer _fp[3];
                char _zp[3];

                AdjTypePack() {
                    _fp[0] = 0;
                    _fp[1] = 0;
                    _fp[2] = 0;
                }
            };

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
        };

        template < class T> class InfoOcf : public T {
        public:
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
                assert(j >= 0 && j < 3);
                return v[j];
            }
            inline typename T::VertexType *cV(const int j) const
            {
                assert(j >= 0 && j < 3);
                return v[j];
            }

            inline typename T::VertexType *&V0(const int j)
            {
                return V(j);
            }
            inline typename T::VertexType *&V1(const int j)
            {
                return V((j + 1) % 3);
            }
            inline typename T::VertexType *&V2(const int j)
            {
                return V((j + 2) % 3);
            }

        private:
            typename T::VertexType *v[3];
        };

        template <class FaceType>
        class Pos
        {
        public:
            typedef typename FaceType::VertexType VertexType;
            typedef Pos<FaceType> PosType;
            typedef typename VertexType::ScalarType ScalarType;
            typename FaceType::FaceType *f;

            int z;
            VertexType *v;

            Pos() : f(0), z(-1), v(0) {}
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

            inline bool operator == (PosType const &p) const
            {
                return (f == p.f && z == p.z && v == p.v);
            }

            inline bool operator != (PosType const &p) const
            {
                return (f != p.f || z != p.z || v != p.v);
            }

            void NextF()
            {
                FaceType *t = f;
                f = t->FFp(z);
                z = t->FFi(z);
            }

            void NextE()
            {
                assert(f->V(z) == v || f->V(f->Next(z)) == v);
                FlipE();
                FlipF();
                assert(f->V(z) == v || f->V(f->Next(z)) == v);
            }

            void FlipE()
            {
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V((z + 0) % f->VN()) == v));
                if (f->V(f->Next(z)) == v) z = f->Next(z);
                else z = f->Prev(z);
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V((z)) == v));
            }

            void FlipF()
            {
                assert(f->FFp(z)->FFp(f->FFi(z)) == f);
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V((z)) == v));
                FaceType *nf = f->FFp(z);
                int nz = f->FFi(z);
                assert(nf->V(nf->Prev(nz)) != v && (nf->V(nf->Next(nz)) == v || nf->V((nz)) == v));
                f = nf;
                z = nz;
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));
            }

            void FlipV()
            {
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));

                if (f->V(f->Next(z)) == v)
                    v = f->V(z);
                else
                    v = f->V(f->Next(z));

                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));
            }

            VertexType *VFlip() const
            {
                assert(f->cV(f->Prev(z)) != v && (f->cV(f->Next(z)) == v || f->cV(z) == v));
                if (f->cV(f->Next(z)) == v)	return f->cV(z);
                else			return f->cV(f->Next(z));
            }

            FaceType *FFlip() const
            {
                FaceType *nf = f->FFp(z);
                return nf;
            }

            void NextB()
            {
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));
                assert(f->FFp(z) == f);

                do
                {
                    NextE();
                }
                while (!IsBorder());

                assert(IsBorder() && (f->V(z) == v || f->V(f->Next(z)) == v));

                FlipV();
                assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));
                assert(f->FFp(z) == f);
            }

            bool IsBorder()const
            {
                return face::IsBorder(*f, z);
            }
        };

        template <class FaceType>
        inline bool IsBorder(FaceType const &f, const int j)
        {
            if (FaceType::HasFFAdjacency())
                return f.cFFp(j) == &f;

            assert(0);
            return true;
        }

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

        template <class FaceType>
        bool PointDistanceBase(
            const FaceType &f,	
            const tessellation::Point3<typename FaceType::ScalarType> & q,
            typename FaceType::ScalarType & dist,
            tessellation::Point3<typename FaceType::ScalarType> & p)
        {
            typedef typename FaceType::ScalarType ScalarType;

            if (f.cN() == Point3<ScalarType>(0, 0, 0))
            {
                Box3<ScalarType> bb;
                f.GetBBox(bb);
                Segment3<ScalarType> degenTri(bb.min, bb.max);
                Point3<ScalarType> closest;
                ScalarType d;
                if (bb.Diag() > 0)
                    tessellation::SegmentPointDistance<ScalarType>(degenTri, q, closest, d);
                else
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

            ScalarType d = SignedDistancePlanePoint(fPlane, q);
            if (d > dist || d < -dist)
                return false;

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
                if (fabs(f.cN()[1]) > fabs(f.cN()[2])) bestAxis = 1;
                else bestAxis = 2;
            }

            ScalarType scaleFactor;

            switch (bestAxis)
            {
            case 0:
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
                    if (dist > bt) { dist = bt; return true; }
                    else return false;
                }
                break;

            case 1:
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
                    if (dist > bt) { dist = bt; return true; }
                    else return false;
                }
                break;

            case 2:
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

                    if (dist > bt) { dist = bt; return true; }
                    else return false;
                }
                break;
            default: assert(0);

            }

            dist = ScalarType(fabs(d));
            return true;
        }

    }

    template <class UserTypes>
    class FaceBase : public
        face::EmptyCore< FaceTypeHolder <UserTypes> >
    {
    };
    
    template < class UserTypes,
        template <typename> class A, template <typename> class B, template <typename> class C,
        template <typename> class D, template <typename> class E,
        template <typename> class F, template <typename> class G>
    class FaceArityMax : public  Arity7<FaceBase<UserTypes>, A, B, C, D, E, F, G>
    {
    public:
        typedef typename  FaceArityMax::ScalarType ScalarType;
        enum
        {
            DELETED = 0x00000001,	
            NOTREAD = 0x00000002,	
            NOTWRITE = 0x00000004,	
            VISITED = 0x00000010,	
            SELECTED = 0x00000020,	
            BORDER0 = 0x00000040,
            BORDER1 = 0x00000080,
            BORDER2 = 0x00000100,
            BORDER012 = BORDER0 | BORDER1 | BORDER2,
            FAUX0 = 0x00040000,
            FAUX1 = 0x00080000,
            FAUX2 = 0x00100000,
            FAUX012 = FAUX0 | FAUX1 | FAUX2,
            USER0 = 0x00200000
        };

        bool IsD() const
        {
            return (this->cFlags() & DELETED) != 0;
        }

        bool IsR() const
        {
            return (this->cFlags() & NOTREAD) == 0;
        }

        bool IsW() const
        {
            return (this->cFlags() & NOTWRITE) == 0;
        }

        bool IsRW() const
        {
            return (this->cFlags() & (NOTREAD | NOTWRITE)) == 0;
        }

        bool IsS() const
        {
            return (this->cFlags() & SELECTED) != 0;
        }

        bool IsV() const
        {
            return (this->cFlags() & VISITED) != 0;
        }

        void SetFlags(int flagp)
        {
            this->Flags() = flagp;
        }

        void ClearFlags()
        {
            this->Flags() = 0;
        }
        
        void SetS()
        {
            this->Flags() |= SELECTED;
        }

        void ClearS()
        {
            this->Flags() &= ~SELECTED;
        }

        void SetV()
        {
            this->Flags() |= VISITED;
        }

        void ClearV()
        {
            this->Flags() &= ~VISITED;
        }

        bool IsB(int i) const
        {
            return (this->cFlags() & (BORDER0 << i)) != 0;
        }

        void SetB(int i)
        {
            this->Flags() |= (BORDER0 << i);
        }

        void ClearB(int i)
        {
            this->Flags() &= (~(BORDER0 << i));
        }

        bool IsF(int i) const
        {
            return (this->cFlags() & (FAUX0 << i)) != 0;
        }
        bool IsAnyF() const
        {
            return (this->cFlags() & (FAUX0 | FAUX1 | FAUX2)) != 0;
        }

        void SetF(int i)
        {
            this->Flags() |= (FAUX0 << i);
        }

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
