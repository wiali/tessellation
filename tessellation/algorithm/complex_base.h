#pragma once

#include <assert.h>
#include <string>
#include <vector>
#include <set>
#include <stack>
#include <algorithm>
#include <map>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <set>
#include <algorithm/space_base.h>
#include <algorithm/simplex_base.h>

namespace tessellation
{

    class MissingComponentException : public std::runtime_error
    {
    public:
        MissingComponentException(const std::string &err) : std::runtime_error(err)
        {
            std::cout << "Missing Component Exception -" << err << "- \n";
        }
        virtual const char *what() const throw ()
        {
            static char buf[128] = "Missing Component";
            return buf;
        }
    };

    // dummy mesh

    struct _Vertex;
    struct _Face;

    struct DummyTypes
    {
        typedef tessellation::Point3<bool> CoordType;
        typedef char ScalarType;
    };

    template <class A>
    struct Use
    {
        template <class T> struct AsVertexType : public T
        {
            typedef A VertexType;
            typedef VertexType *VertexPointer;
        };
        template <class T> struct AsFaceType : public T
        {
            typedef A FaceType;
            typedef FaceType *FacePointer;
        };
    };

    template < template <typename> class A = DefaultDeriver, template <typename> class B = DefaultDeriver,
             template <typename> class C = DefaultDeriver >
    class UsedTypes : public Arity5 < DummyTypes,
        Use<  Vertex	<UsedTypes< A, B, C> > > :: template AsVertexType,
              Use<  Face		<UsedTypes< A, B, C> > > :: template AsFaceType,
                        A, B, C
                        >
    {
    };

    struct _UsedTypes : public UsedTypes <
            Use<_Vertex>::AsVertexType,
            Use<_Face  >::AsFaceType
            > {};

    struct _Vertex : public  Vertex<_UsedTypes> {};
    struct _Face : public  Face<_UsedTypes> {};


    class SimpleTempDataBase
    {
    public:
        virtual ~SimpleTempDataBase() {}
        SimpleTempDataBase() {}
        virtual void Resize(size_t sz) = 0;
        virtual void Reorder(std::vector<size_t> &newVertIndex) = 0;
    };

    template <class TYPE>
    class VectorNBW : public std::vector<TYPE> {};


    template <class STL_CONT, class ATTR_TYPE>
    class SimpleTempData: public SimpleTempDataBase
    {

    public:
        typedef SimpleTempData<STL_CONT, ATTR_TYPE> SimpTempDataType;
        typedef ATTR_TYPE AttrType;

        STL_CONT &c;
        VectorNBW<ATTR_TYPE> data;
        int padding;

        SimpleTempData(STL_CONT  &_c): c(_c), padding(0)
        {
            data.reserve(c.capacity());
            data.resize(c.size());
        };
        SimpleTempData(STL_CONT  &_c, const ATTR_TYPE &val): c(_c)
        {
            data.reserve(c.capacity());
            data.resize(c.size());
            Init(val);
        };

        ~SimpleTempData()
        {
            data.clear();
        }

        void Init(const ATTR_TYPE &val)
        {
            std::fill(data.begin(), data.end(), val);
        }
        // access to data
        ATTR_TYPE &operator[](const typename STL_CONT::value_type &v)
        {
            return data[&v - &*c.begin()];
        }
        ATTR_TYPE &operator[](const typename STL_CONT::value_type *v)
        {
            return data[v - &*c.begin()];
        }
        ATTR_TYPE &operator[](const typename STL_CONT::iterator &cont)
        {
            return data[&(*cont) - &*c.begin()];
        }
        ATTR_TYPE &operator[](size_t i)
        {
            return data[i];
        }


        void Resize(size_t sz)
        {
            data.resize(sz);
        }

        void Reorder(std::vector<size_t> &newVertIndex)
        {
            for(unsigned int i = 0 ; i < data.size(); ++i)
            {
                if( newVertIndex[i] != (std::numeric_limits<size_t>::max)())
                    data[newVertIndex[i]] = data[i];
            }
        }

    };

    template <class ATTR_TYPE>
    class Attribute: public SimpleTempDataBase
    {
    public:
        typedef ATTR_TYPE AttrType;
        AttrType *attribute;
        Attribute()
        {
            attribute = new ATTR_TYPE();
        }
        ~Attribute()
        {
            delete attribute;
        }

        void Resize(size_t  )
        {
            assert(0);
        }
        void Reorder(std::vector<size_t> &  )
        {
            assert(0);
        }

    };

    class PointerToAttribute
    {
    public:
        SimpleTempDataBase *_handle;		// pointer to the SimpleTempData that stores the attribute
        std::string _name;					// name of the attribute
        int _sizeof;						// size of the attribute type (used only with VMI loading)
        int _padding;						// padding 	(used only with VMI loading)

        int n_attr;							// unique ID of the attribute

        void Resize(size_t sz)
        {
            ((SimpleTempDataBase *)_handle)->Resize(sz);
        }
        void Reorder(std::vector<size_t> &newVertIndex)
        {
            ((SimpleTempDataBase *)_handle)->Reorder(newVertIndex);
        }
        bool operator<(const  PointerToAttribute    b) const
        {
            return(_name.empty() && b._name.empty()) ? (_handle < b._handle) : ( _name < b._name);
        }
    };


    namespace tri
    {
        /** \addtogroup trimesh */
        /*@{*/


        /* MeshTypeHolder is a class which is used to define the types in the mesh
        */

        template <class TYPESPOOL>
        struct BaseMeshTypeHolder
        {

            typedef bool ScalarType;
            typedef std::vector< typename TYPESPOOL::VertexType  >	CONTV;
            typedef std::vector< typename TYPESPOOL::FaceType >		CONTF;

            typedef CONTV									VertContainer;
            typedef _Vertex 								VertexType;
            typedef typename TYPESPOOL::VertexPointer		VertexPointer;
            typedef const typename TYPESPOOL::VertexPointer ConstVertexPointer;
            typedef bool									CoordType;
            typedef typename CONTV::iterator				VertexIterator;
            typedef typename CONTV::const_iterator	ConstVertexIterator;


            typedef CONTF														FaceContainer;
            typedef typename CONTF::value_type					FaceType;
            typedef typename CONTF::const_iterator				ConstFaceIterator;
            typedef typename CONTF::iterator					FaceIterator;
            typedef typename TYPESPOOL::FacePointer	FacePointer;
            typedef const typename TYPESPOOL::FacePointer		ConstFacePointer;

        };
        
        template <class T, typename CONT, class TRAIT >
        struct MeshTypeHolder: public T {};

        template <class T, typename CONT>
        struct MeshTypeHolder<T, CONT, AllTypes::AVertexType>: public T
        {
            typedef CONT VertContainer;
            typedef typename VertContainer::value_type VertexType;
            typedef VertexType *VertexPointer;
            typedef const VertexType *ConstVertexPointer;
            typedef typename VertexType::ScalarType ScalarType;
            typedef typename VertexType::CoordType CoordType;
            typedef typename VertContainer::iterator VertexIterator;
            typedef typename VertContainer::const_iterator ConstVertexIterator;
        };

        template <typename T, class CONT>
        struct MeshTypeHolder< T, CONT, AllTypes::AFaceType> :public T {
            typedef CONT FaceContainer;
            typedef typename FaceContainer::value_type FaceType;
            typedef typename FaceContainer::const_iterator ConstFaceIterator;
            typedef typename FaceContainer::iterator FaceIterator;
            typedef FaceType * FacePointer;
            typedef const FaceType * ConstFacePointer;
        };

        template <typename T, typename CONT> struct Der : public MeshTypeHolder<T, CONT, typename CONT::value_type::IAm> {};
        struct DummyContainer
        {
            struct value_type
            {
                typedef int IAm;
            };
        };


        template < class Container0 = DummyContainer, class Container1 = DummyContainer, class Container2 = DummyContainer, class Container3 = DummyContainer >
        class TriMesh
            : public  MArity4<BaseMeshTypeHolder<typename Container0::value_type::TypesPool>, Container0, Der, Container1, Der, Container2, Der, Container3, Der>
        {
        public:

            typedef typename TriMesh::ScalarType		ScalarType;
            typedef typename TriMesh::VertContainer VertContainer;
            typedef typename TriMesh::FaceContainer FaceContainer;

            // types for vertex
            typedef typename TriMesh::VertexType						VertexType;
            typedef typename TriMesh::VertexPointer					VertexPointer;
            typedef typename TriMesh::ConstVertexPointer		ConstVertexPointer;
            typedef typename TriMesh::CoordType							CoordType;
            typedef typename TriMesh::VertexIterator				VertexIterator;
            typedef typename TriMesh::ConstVertexIterator		ConstVertexIterator;


            //types for face
            typedef typename TriMesh::FaceType							FaceType;
            typedef typename TriMesh::ConstFaceIterator			ConstFaceIterator;
            typedef typename TriMesh::FaceIterator					FaceIterator;
            typedef typename TriMesh::FacePointer						FacePointer;
            typedef typename TriMesh::ConstFacePointer			ConstFacePointer;

            typedef tessellation::PointerToAttribute PointerToAttribute;

            typedef TriMesh<Container0, Container1> MeshType;

            typedef Box3<ScalarType> BoxType;

            /// Container of vertices, usually a vector.
            VertContainer vert;
            /// Current number of vertices; this member is for internal use only. You should always use the VN() member
            int vn;
            /// Current number of vertices
            inline int VN() const
            {
                return vn;
            }
            
            /// Container of faces, usually a vector.
            FaceContainer face;
            /// Current number of faces; this member is for internal use only. You should always use the FN() member
            int fn;
            /// Current number of faces
            inline int FN() const
            {
                return fn;
            }

            /// Nomi di textures
            //
            std::vector<std::string> textures;
            
            int attrn;	// total numer of attribute created
            
        /// Bounding box of the mesh
            Box3<typename TriMesh::VertexType::CoordType::ScalarType> bbox;

            std::set< PointerToAttribute > vert_attr;
            std::set< PointerToAttribute > face_attr;
            std::set< PointerToAttribute > mesh_attr;
            
            template <class ATTR_TYPE, class CONT>
            class AttributeHandle
            {
            public:
                AttributeHandle()
                {
                    _handle = (SimpleTempData<CONT, ATTR_TYPE> *)NULL;
                }
                AttributeHandle( void *ah, const int &n): _handle ( (SimpleTempData<CONT, ATTR_TYPE> *)ah ), n_attr(n) {}
                AttributeHandle operator = ( const PointerToAttribute &pva)
                {
                    _handle = (SimpleTempData<CONT, ATTR_TYPE> *)pva._handle;
                    n_attr = pva.n_attr;
                    return (*this);
                }

                //pointer to the SimpleTempData that stores the attribute
                SimpleTempData<CONT, ATTR_TYPE> *_handle;

                // its attribute number
                int n_attr;

                // access function
                template <class RefType>
                ATTR_TYPE &operator [](const RefType   &i)
                {
                    return (*_handle)[i];
                }
                virtual void resize(size_t /*size*/) { };
            };

            template <class ATTR_TYPE>
            class PerVertexAttributeHandle: public AttributeHandle<ATTR_TYPE, VertContainer>
            {
            public:
                PerVertexAttributeHandle(): AttributeHandle<ATTR_TYPE, VertContainer>() {}
                PerVertexAttributeHandle( void *ah, const int &n): AttributeHandle<ATTR_TYPE, VertContainer>(ah, n) {}
            };

            template <class ATTR_TYPE>
            class PerFaceAttributeHandle : public AttributeHandle<ATTR_TYPE, FaceContainer>
            {
            public:
                PerFaceAttributeHandle() : AttributeHandle<ATTR_TYPE, FaceContainer>() {}
                PerFaceAttributeHandle(void *ah, const int &n) : AttributeHandle<ATTR_TYPE, FaceContainer>(ah, n) {}
            };


            template <class ATTR_TYPE>
            class PerMeshAttributeHandle
            {
            public:
                PerMeshAttributeHandle()
                {
                    _handle = NULL;
                }
                PerMeshAttributeHandle(void *ah, const int &n): _handle ( (Attribute<ATTR_TYPE> *)ah ), n_attr(n) {}
                PerMeshAttributeHandle operator = ( const PerMeshAttributeHandle &pva)
                {
                    _handle = (Attribute<ATTR_TYPE> *)pva._handle;
                    n_attr = pva.n_attr;
                    return (*this);
                }

                Attribute<ATTR_TYPE> *_handle;
                int n_attr;
                ATTR_TYPE &operator ()()
                {
                    return *((Attribute<ATTR_TYPE> *)_handle)->attribute;
                }
            };

        private:
            /// The per-mesh color. Not very useful and meaningful...
            Color4b c;
        public:

            inline const Color4b &C() const
            {
                return c;
            }
            inline       Color4b &C()
            {
                return c;
            }
            inline       Color4b cC() const
            {
                return c;
            }

            /// Default constructor
            TriMesh()
            {
                Clear();
            }

            /// destructor
            ~TriMesh()
            {
                Clear();
            }

            /// Function to destroy the mesh
            void Clear()
            {
                for(FaceIterator fi = face.begin(); fi != face.end(); ++fi)
                    (*fi).Dealloc();
                vert.clear();
                face.clear();
                //    textures.clear();
                //    normalmaps.clear();
                vn = 0;
                fn = 0;
                imark = 0;
                C() = Color4b::Gray;
            }
            
            bool IsEmpty() const
            {
                return vert.empty() && face.empty();
            }

            int &SimplexNumber()
            {
                return fn;
            }
            int &VertexNumber()
            {
                return vn;
            }

            /// The incremental mark
            int imark;

        private:
            // TriMesh cannot be copied. Use Append (see vcg/complex/append.h)
            TriMesh operator =(const TriMesh &  /*m*/)
            {
                assert(0);
                return TriMesh();
            }
            TriMesh(const TriMesh & ) {}

        };	// end class Mesh
        
        template < class VertexType> bool VertexVectorHasVFAdjacency(const std::vector<VertexType> &)
        {
            return VertexType::HasVFAdjacency();
        }
        template < class VertexType> bool VertexVectorHasVEAdjacency(const std::vector<VertexType> &)
        {
            return VertexType::HasVEAdjacency();
        }

        template < class FaceType  > bool   FaceVectorHasVFAdjacency(const std::vector<FaceType  > &)
        {
            return FaceType::HasVFAdjacency();
        }

        template < class TriMeshType> bool HasPerVertexVFAdjacency     (const TriMeshType &m)
        {
            return tri::VertexVectorHasVFAdjacency(m.vert);
        }
        template < class TriMeshType> bool HasPerVertexVEAdjacency     (const TriMeshType &m)
        {
            return tri::VertexVectorHasVEAdjacency(m.vert);
        }

        template < class TriMeshType> bool   HasPerFaceVFAdjacency     (const TriMeshType &m)
        {
            return tri::FaceVectorHasVFAdjacency  (m.face);
        }
        
        template < class VertexType> bool VertexVectorHasPerVertexFlags       (const std::vector<VertexType> &)
        {
            return VertexType::HasFlags       ();
        }

        template < class TriMeshType> bool HasPerVertexFlags       (const TriMeshType &m)
        {
            return tri::VertexVectorHasPerVertexFlags       (m.vert);
        }

        template < class FaceType>    bool FaceVectorHasPerFaceFlags  (const std::vector<FaceType> &)
        {
            return FaceType::HasFlags  ();
        }
        template < class FaceType>    bool FaceVectorHasPerFaceNormal (const std::vector<FaceType> &)
        {
            return FaceType::HasNormal ();
        }
        template < class FaceType>    bool FaceVectorHasPerFaceColor  (const std::vector<FaceType> &)
        {
            return FaceType::HasColor  ();
        }
        template < class FaceType>    bool FaceVectorHasPerFaceMark   (const std::vector<FaceType> &)
        {
            return FaceType::HasMark   ();
        }
        template < class FaceType>    bool FaceVectorHasFFAdjacency   (const std::vector<FaceType> &)
        {
            return FaceType::HasFFAdjacency();
        }
        template < class FaceType>    bool FaceVectorHasFEAdjacency   (const std::vector<FaceType> &)
        {
            return FaceType::HasFEAdjacency();
        }
        template < class FaceType>    bool FaceVectorHasFVAdjacency   (const std::vector<FaceType> &)
        {
            return FaceType::HasFVAdjacency();
        }

        template < class TriMeshType> bool HasPerFaceNormal(const TriMeshType &m) 
        { 
            return tri::FaceVectorHasPerFaceNormal(m.face); 
        }

        template < class TriMeshType> bool HasPerVertexNormal(const TriMeshType &m)
        { 
            return tri::VertexVectorHasPerVertexNormal(m.vert); 
        }

        template < class VertexType >
        bool VertexVectorHasPerVertexNormal(const vertex::vector_ocf<VertexType> &fv)
        {
            if (VertexType::HasNormalOcf()) return fv.IsNormalEnabled();
            else return VertexType::HasNormal();
        }

        template < class TriMeshType> bool HasPerFaceFlags       (const TriMeshType &m)
        {
            return tri::FaceVectorHasPerFaceFlags       (m.face);
        }
        template < class TriMeshType> bool HasPerFaceColor       (const TriMeshType &m)
        {
            return tri::FaceVectorHasPerFaceColor       (m.face);
        }
        template < class TriMeshType> bool HasFFAdjacency   (const TriMeshType &m)
        {
            return tri::FaceVectorHasFFAdjacency   (m.face);
        }
        template < class TriMeshType> bool HasFEAdjacency   (const TriMeshType &m)
        {
            return tri::FaceVectorHasFEAdjacency   (m.face);
        }
        template < class TriMeshType> bool HasFVAdjacency   (const TriMeshType &m)
        {
            return tri::FaceVectorHasFVAdjacency   (m.face);
        }

        template < class TriMeshType> bool HasVFAdjacency   (const TriMeshType &m)
        {
            return tri::FaceVectorHasVFAdjacency(m.face) 
                && tri::VertexVectorHasVFAdjacency(m.vert);
        }

        template < class  CType0, class CType1>
        bool HasVHAdjacency (const TriMesh < CType0, CType1> & /*m*/)
        {
            return TriMesh < CType0 , CType1>::VertContainer::value_type::HasVHAdjacency();
        }
        
        template < class TriMeshType> bool HasPerFaceMark(const TriMeshType &m) { return tri::FaceVectorHasPerFaceMark(m.face); }

        template <class MeshType> void RequirePerVertexFlags(MeshType &m) 
        {
            if (!tri::HasPerVertexFlags(m)) throw tessellation::MissingComponentException("PerVertexFlags       "); 
        }

        template <class MeshType> void RequireVFAdjacency(MeshType &m)
        {
            if (!tri::HasVFAdjacency(m)) throw tessellation::MissingComponentException("VFAdjacency");
        }
        template <class MeshType> void RequireVEAdjacency(MeshType &m)
        {
            if (!tri::HasVEAdjacency(m)) throw tessellation::MissingComponentException("VEAdjacency");
        }
        template <class MeshType> void RequireFFAdjacency(MeshType &m)
        {
            if (!tri::HasFFAdjacency(m)) throw tessellation::MissingComponentException("FFAdjacency");
        }
        template <class MeshType> void RequireEEAdjacency(MeshType &m)
        {
            if (!tri::HasEEAdjacency(m)) throw tessellation::MissingComponentException("EEAdjacency");
        }
        template <class MeshType> void RequireFEAdjacency(MeshType &m)
        {
            if (!tri::HasFEAdjacency(m)) throw tessellation::MissingComponentException("FEAdjacency");
        }
        template <class MeshType> void RequireFHAdjacency(MeshType &m)
        {
            if (!tri::HasFHAdjacency(m)) throw tessellation::MissingComponentException("FHAdjacency");
        }

        template <class MeshType> void RequirePerFaceFlags(MeshType &m)
        {
            if (!tri::HasPerFaceFlags(m)) throw tessellation::MissingComponentException("PerFaceFlags       ");
        }

        template <class MeshType> void RequirePerFaceNormal(MeshType &m) 
        { 
            if (!tri::HasPerFaceNormal(m)) throw tessellation::MissingComponentException("PerFaceNormal      "); 
        }

        template <class MeshType> void RequirePerVertexNormal(MeshType &m) 
        { 
            if (!tri::HasPerVertexNormal(m)) throw tessellation::MissingComponentException("PerVertexNormal      ");
        }
        
        template <class MeshType> void RequirePerFaceMark(MeshType &m) { if (!tri::HasPerFaceMark(m)) throw tessellation::MissingComponentException("PerFaceMark        "); }


        template<class MeshType>
        size_t Index(MeshType &m, const typename MeshType::VertexType *vp)
        {
            return vp - &*m.vert.begin();
        }
        template<class MeshType>
        size_t Index(MeshType &m, const typename MeshType::FaceType *fp)
        {
            return fp - &*m.face.begin();
        }

        template <class MeshType, class ATTR_CONT>
        void ReorderAttribute(ATTR_CONT &c, std::vector<size_t> & newVertIndex, MeshType & /* m */) {
            typename std::set<typename MeshType::PointerToAttribute>::iterator ai;
            for (ai = c.begin(); ai != c.end(); ++ai)
                ((typename MeshType::PointerToAttribute)(*ai)).Reorder(newVertIndex);
        }

        template <class MeshType, class ATTR_CONT>
        void ResizeAttribute(ATTR_CONT &c, size_t sz, MeshType &/*m*/) {
            typename std::set<typename MeshType::PointerToAttribute>::iterator ai;
            for (ai = c.begin(); ai != c.end(); ++ai)
                ((typename MeshType::PointerToAttribute)(*ai)).Resize(sz);
        }

        
        /*!
        \brief  Class to safely add and delete elements in a mesh.

        Adding elements to a mesh, like faces and vertices can involve the reallocation of the vectors of the involved elements.
        This class provide the only safe methods to add elements.
        It also provide an accessory class tessellation::tri::PointerUpdater for updating pointers to mesh elements that are kept by the user.
        */
        template <class MeshType>
        class Allocator
        {

        public:
            typedef typename MeshType::VertexType     VertexType;
            typedef typename MeshType::VertexPointer  VertexPointer;
            typedef typename MeshType::VertexIterator VertexIterator;
            typedef typename MeshType::VertContainer VertContainer;

            typedef typename MeshType::FaceType       FaceType;
            typedef typename MeshType::FacePointer    FacePointer;
            typedef typename MeshType::FaceIterator   FaceIterator;
            typedef typename MeshType::FaceContainer FaceContainer;

            typedef typename MeshType::CoordType     CoordType;


            typedef typename MeshType::PointerToAttribute PointerToAttribute;
            typedef typename std::set<PointerToAttribute>::iterator AttrIterator;
            typedef typename std::set<PointerToAttribute>::const_iterator AttrConstIterator;
            typedef typename std::set<PointerToAttribute >::iterator PAIte;

            template<class SimplexPointerType>
            class PointerUpdater
            {
            public:
                PointerUpdater(void) : newBase(0), oldBase(0), newEnd(0), oldEnd(0), preventUpdateFlag(false)
                {
                    ;
                }
                void Clear()
                {
                    newBase = oldBase = newEnd = oldEnd = 0;
                    remap.clear();
                }
                /*! \brief Update a pointer to an element of a mesh after a reallocation

                The updating is correctly done only if this PointerUpdater have been passed to the corresponing allocation call. \sa \ref allocation
                */
                void Update(SimplexPointerType &vp)
                {
                    //if(vp>=newBase && vp<newEnd) return;
                    if (vp < oldBase || vp > oldEnd) return;
                    assert(vp >= oldBase);
                    assert(vp < oldEnd);
                    vp = newBase + (vp - oldBase);
                    if (!remap.empty())
                        vp = newBase + remap[vp - newBase];
                }
                /*!
                \brief return true if the allocation operation that initialized this PointerUpdater has caused a reallocation
                */
                bool NeedUpdate()
                {
                    if ((oldBase && newBase != oldBase && !preventUpdateFlag) || !remap.empty()) return true;
                    else return false;
                }

                SimplexPointerType newBase;
                SimplexPointerType oldBase;
                SimplexPointerType newEnd;
                SimplexPointerType oldEnd;
                std::vector<size_t> remap; // this vector keep the new position of an element. Uninitialized elements have max_int value to denote an element that has not to be remapped.

                bool preventUpdateFlag; /// when true no update is considered necessary.
            };

            /* +++++++++++++++ Add Vertices ++++++++++++++++ */

            /** \brief Add n vertices to the mesh.
            Function to add n vertices to the mesh.
            The elements are added always to the end of the vector.
            No attempt of reusing previously deleted element is done.
            \sa PointerUpdater
            \param m the mesh to be modified
            \param n the number of elements to be added
            \param pu  a PointerUpdater initialized so that it can be used to update pointers to vertices that could have become invalid after this adding.
            \retval the iterator to the first element added.
            */
            static VertexIterator AddVertices(MeshType &m, size_t n, PointerUpdater<VertexPointer> &pu)
            {
                VertexIterator last;
                if (n == 0) return m.vert.end();
                pu.Clear();
                if (m.vert.empty()) pu.oldBase = 0;  // if the vector is empty we cannot find the last valid element
                else
                {
                    pu.oldBase = &*m.vert.begin();
                    pu.oldEnd = &m.vert.back() + 1;
                }

                m.vert.resize(m.vert.size() + n);
                m.vn += int(n);

                typename std::set<PointerToAttribute>::iterator ai;
                for (ai = m.vert_attr.begin(); ai != m.vert_attr.end(); ++ai)
                    ((PointerToAttribute)(*ai)).Resize(m.vert.size());

                pu.newBase = &*m.vert.begin();
                pu.newEnd = &m.vert.back() + 1;
                if (pu.NeedUpdate())
                {
                    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
                        if (!(*fi).IsD())
                            for (int i = 0; i < (*fi).VN(); ++i)
                                if ((*fi).cV(i) != 0) pu.Update((*fi).V(i));

                }
                size_t siz = (size_t)(m.vert.size() - n);

                last = m.vert.begin();
                advance(last, siz);

                return last;// deve restituire l'iteratore alla prima faccia aggiunta;
            }

            /** \brief Wrapper to AddVertices(); no PointerUpdater
            */
            static VertexIterator AddVertices(MeshType &m, size_t n)
            {
                PointerUpdater<VertexPointer> pu;
                return AddVertices(m, n, pu);
            }
            
            /* +++++++++++++++ Add Faces ++++++++++++++++ */

            /** Function to add a face to the mesh and initializing it with the three given VertexPointers
            */
            static FaceIterator AddFace(MeshType &m, VertexPointer v0, VertexPointer v1, VertexPointer v2)
            {
                assert(m.vert.size() > 0);
                assert((v0 != v1) && (v1 != v2) && (v0 != v2));
                assert(v0 >= &m.vert.front() && v0 <= &m.vert.back());
                assert(v1 >= &m.vert.front() && v1 <= &m.vert.back());
                assert(v2 >= &m.vert.front() && v2 <= &m.vert.back());
                PointerUpdater<FacePointer> pu;
                FaceIterator fi = AddFaces(m, 1, pu);
                fi->Alloc(3);
                fi->V(0) = v0;
                fi->V(1) = v1;
                fi->V(2) = v2;
                return fi;
            }

            /** Function to add a face to the mesh and initializing it with three indexes
            */
            static FaceIterator AddFace(MeshType &m, size_t v0, size_t v1, size_t v2)
            {
                assert((v0 != v1) && (v1 != v2) && (v0 != v2));
                assert(v0 >= 0 && v0 <= m.vert.size());
                assert(v1 >= 0 && v1 <= m.vert.size());
                assert(v2 >= 0 && v2 <= m.vert.size());
                return AddFace(m, &(m.vert[v0]), &(m.vert[v1]), &(m.vert[v2]));
            }


            /** \brief Function to add n faces to the mesh.
            First wrapper, with no parameters
            */
            static FaceIterator AddFaces(MeshType &m, size_t n)
            {
                PointerUpdater<FacePointer> pu;
                return AddFaces(m, n, pu);
            }

            /** \brief Function to add n faces to the mesh.
            Second Wrapper, with a vector of face pointer to be updated.
            */
            static FaceIterator AddFaces(MeshType &m, size_t n, std::vector<FacePointer *> &local_vec)
            {
                PointerUpdater<FacePointer> pu;
                FaceIterator f_ret = AddFaces(m, n, pu);

                typename std::vector<FacePointer *>::iterator fi;
                for (fi = local_vec.begin(); fi != local_vec.end(); ++fi)
                    pu.Update(**fi);
                return f_ret;
            }

            static FaceIterator AddFaces(MeshType &m, size_t n, PointerUpdater<FacePointer> &pu)
            {
                pu.Clear();
                if (n == 0) return m.face.end();
                if (!m.face.empty()) // if the vector is empty we cannot find the last valid element
                {
                    pu.oldBase = &*m.face.begin();
                    pu.oldEnd = &m.face.back() + 1;
                }
                // The actual resize
                m.face.resize(m.face.size() + n);
                m.fn += int(n);

                size_t siz = (size_t)(m.face.size() - n);
                FaceIterator firstNewFace = m.face.begin();
                advance(firstNewFace, siz);

                typename std::set<PointerToAttribute>::iterator ai;
                for (ai = m.face_attr.begin(); ai != m.face_attr.end(); ++ai)
                    ((PointerToAttribute)(*ai)).Resize(m.face.size());

                pu.newBase = &*m.face.begin();
                pu.newEnd = &m.face.back() + 1;

                if (pu.NeedUpdate())
                {
                    if (HasFFAdjacency(m))
                    {
                        // cycle on all the faces except the new ones
                        for (FaceIterator fi = m.face.begin(); fi != firstNewFace; ++fi)
                            if (!(*fi).IsD())
                                for (int i = 0; i < (*fi).VN(); ++i)
                                    if ((*fi).cFFp(i) != 0) pu.Update((*fi).FFp(i));
                    }

                    if (HasPerVertexVFAdjacency(m) && HasPerFaceVFAdjacency(m))
                    {
                        // cycle on all the faces except the new ones
                        for (FaceIterator fi = m.face.begin(); fi != firstNewFace; ++fi)
                            if (!(*fi).IsD())
                                for (int i = 0; i < (*fi).VN(); ++i)
                                    if ((*fi).cVFp(i) != 0) pu.Update((*fi).VFp(i));

                        for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                            if (!(*vi).IsD() && (*vi).cVFp() != 0)
                                pu.Update((*vi).VFp());
                    }
                }
                return firstNewFace;
            }


        public:
            /*! \brief Add a Per-Vertex Attribute of the given ATTR_TYPE with the given name.

            No attribute with that name must exists (even of different type)
            */
            template <class ATTR_TYPE>
            static
            typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
            AddPerVertexAttribute(MeshType &m, std::string name)
            {
                PAIte i;
                PointerToAttribute h;
                h._name = name;
                if (!name.empty())
                {
                    i = m.vert_attr.find(h);
                    assert(i == m.vert_attr.end());// an attribute with this name exists
                }

                h._sizeof = sizeof(ATTR_TYPE);
                h._padding = 0;
                h._handle = new SimpleTempData<VertContainer, ATTR_TYPE>(m.vert);
                m.attrn++;
                h.n_attr = m.attrn;
                std::pair < AttrIterator, bool> res = m.vert_attr.insert(h);
                return typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>(res.first->_handle, res.first->n_attr);
            }

            template <class ATTR_TYPE>
            static typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
            AddPerVertexAttribute(MeshType &m)
            {
                return AddPerVertexAttribute<ATTR_TYPE>(m, std::string(""));
            }

            /*! \brief If  the per-vertex attribute exists, delete it.
            */
            template <class ATTR_TYPE>
            static
            void
            DeletePerVertexAttribute(MeshType &m, typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE> &h)
            {
                typename std::set<PointerToAttribute > ::iterator i;
                for (i = m.vert_attr.begin(); i != m.vert_attr.end(); ++i)
                    if ((*i)._handle == h._handle)
                    {
                        delete ((SimpleTempData<VertContainer, ATTR_TYPE> *)(*i)._handle);
                        m.vert_attr.erase(i);
                        return;
                    }
            }

            /// Per Face Attributes
            template <class ATTR_TYPE>
            static
            typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>
            AddPerFaceAttribute(MeshType &m, std::string name)
            {
                PAIte i;
                PointerToAttribute h;
                h._name = name;
                if (!name.empty())
                {
                    i = m.face_attr.find(h);
                    assert(i == m.face_attr.end());// an attribute with this name exists
                }

                h._sizeof = sizeof(ATTR_TYPE);
                h._padding = 0;
                h._handle = new SimpleTempData<FaceContainer, ATTR_TYPE>(m.face);
                m.attrn++;
                h.n_attr = m.attrn;
                std::pair < AttrIterator, bool> res = m.face_attr.insert(h);
                return typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>(res.first->_handle, res.first->n_attr);
            }

            /*! \brief If  the per-face attribute exists, delete it.
            */
            template <class ATTR_TYPE>
            static void DeletePerFaceAttribute(MeshType &m, typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE> &h)
            {
                typename std::set<PointerToAttribute > ::iterator i;
                for (i = m.face_attr.begin(); i != m.face_attr.end(); ++i)
                    if ((*i)._handle == h._handle)
                    {
                        delete ((SimpleTempData<FaceContainer, ATTR_TYPE> *)(*i)._handle);
                        m.face_attr.erase(i);
                        return;
                    }

            }

            static void CompactFaceVector(MeshType &m)
            {
                PointerUpdater<FacePointer>  pu;
                CompactFaceVector(m, pu);
            }

            /*!
            \brief Compact face vector by removing deleted elements.

            Deleted elements are put to the end of the vector and the vector is resized.
            Order between elements is preserved, but not their position (hence the PointerUpdater)
            Immediately after calling this function the \c IsD() test during the scanning a vector, is no more necessary.
            \warning It should not be called when some TemporaryData is active (but works correctly if attributes are present)
            */
            static void CompactFaceVector(MeshType &m, PointerUpdater<FacePointer> &pu)
            {
                // If already compacted fast return please!
                if (m.fn == (int)m.face.size()) return;

                // newFaceIndex [ <old_face_position> ] gives you the new position of the face in the vector;
                pu.remap.resize(m.face.size(), std::numeric_limits<size_t>::max());

                size_t pos = 0;
                for (size_t i = 0; i<m.face.size(); ++i)
                {
                    if (!m.face[i].IsD())
                    {
                        if (pos != i)
                        {
                            m.face[pos].ImportData(m.face[i]);
                            if (FaceType::HasPolyInfo())
                            {
                                m.face[pos].Dealloc();
                                m.face[pos].Alloc(m.face[i].VN());
                            }
                            for (int j = 0; j<m.face[i].VN(); ++j)
                                m.face[pos].V(j) = m.face[i].V(j);

                            if (HasVFAdjacency(m))
                                for (int j = 0; j<m.face[i].VN(); ++j)
                                {
                                    if (m.face[i].IsVFInitialized(j)) {
                                        m.face[pos].VFp(j) = m.face[i].cVFp(j);
                                        m.face[pos].VFi(j) = m.face[i].cVFi(j);
                                    }
                                    else m.face[pos].VFClear(j);
                                }
                            if (HasFFAdjacency(m))
                                for (int j = 0; j<m.face[i].VN(); ++j)
                                {
                                    m.face[pos].FFp(j) = m.face[i].cFFp(j);
                                    m.face[pos].FFi(j) = m.face[i].cFFi(j);
                                }
                        }
                        pu.remap[i] = pos;
                        ++pos;
                    }
                }
                assert((int)pos == m.fn);

                // reorder the optional atttributes in m.face_attr to reflect the changes
                ReorderAttribute(m.face_attr, pu.remap, m);

                FacePointer fbase = &m.face[0];

                // Loop on the vertices to correct VF relation
                if (HasVFAdjacency(m))
                {
                    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                        if (!(*vi).IsD())
                        {
                            if ((*vi).IsVFInitialized() && (*vi).VFp() != 0)
                            {
                                size_t oldIndex = (*vi).cVFp() - fbase;
                                assert(fbase <= (*vi).cVFp() && oldIndex < pu.remap.size());
                                (*vi).VFp() = fbase + pu.remap[oldIndex];
                            }
                        }
                }

                // Loop on the faces to correct VF and FF relations
                pu.oldBase = &m.face[0];
                pu.oldEnd = &m.face.back() + 1;
                for (size_t i = m.fn; i<m.face.size(); ++i)
                    m.face[i].Dealloc();
                m.face.resize(m.fn);
                pu.newBase = (m.face.empty()) ? 0 : &m.face[0];
                pu.newEnd = (m.face.empty()) ? 0 : &m.face.back() + 1;


                // resize the optional atttributes in m.face_attr to reflect the changes
                ResizeAttribute(m.face_attr, m.fn, m);

                // now we update the various (not null) face pointers (inside VF and FF relations)
                for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
                    if (!(*fi).IsD())
                    {
                        if (HasVFAdjacency(m))
                            for (int i = 0; i<(*fi).VN(); ++i)
                                if ((*fi).IsVFInitialized(i) && (*fi).VFp(i) != 0)
                                {
                                    size_t oldIndex = (*fi).VFp(i) - fbase;
                                    assert(fbase <= (*fi).VFp(i) && oldIndex < pu.remap.size());
                                    (*fi).VFp(i) = fbase + pu.remap[oldIndex];
                                }
                        if (HasFFAdjacency(m))
                            for (int i = 0; i<(*fi).VN(); ++i)
                                if ((*fi).cFFp(i) != 0)
                                {
                                    size_t oldIndex = (*fi).FFp(i) - fbase;
                                    assert(fbase <= (*fi).FFp(i) && oldIndex < pu.remap.size());
                                    (*fi).FFp(i) = fbase + pu.remap[oldIndex];
                                }
                    }
            }

            /*! \brief Wrapper without the PointerUpdater. */
            static void CompactVertexVector(MeshType &m) {
                PointerUpdater<VertexPointer>  pu;
                CompactVertexVector(m, pu);
            }

            /*!
            \brief Compact vector of vertices removing deleted elements.
            Deleted elements are put to the end of the vector and the vector is resized. Order between elements is preserved but not their position (hence the PointerUpdater)
            After calling this function the \c IsD() test in the scanning a vector, is no more necessary.

            \warning It should not be called when TemporaryData is active (but works correctly if attributes are present)
            */
            static void CompactVertexVector(MeshType &m, PointerUpdater<VertexPointer> &pu)
            {
                // If already compacted fast return please!
                if (m.vn == (int)m.vert.size()) return;

                // newVertIndex [ <old_vert_position> ] gives you the new position of the vertex in the vector;
                pu.remap.resize(m.vert.size(), std::numeric_limits<size_t>::max());

                size_t pos = 0;
                size_t i = 0;

                for (i = 0; i<m.vert.size(); ++i)
                {
                    if (!m.vert[i].IsD())
                    {
                        pu.remap[i] = pos;
                        ++pos;
                    }
                }
                assert((int)pos == m.vn);

                PermutateVertexVector(m, pu);

            }

            /*
            Function to rearrange the vertex vector according to a given index permutation
            the permutation is vector such that after calling this function

            m.vert[ newVertIndex[i] ] = m.vert[i];

            e.g. newVertIndex[i] is the new index of the vertex i

            */
            static void PermutateVertexVector(MeshType &m, PointerUpdater<VertexPointer> &pu)
            {
                if (m.vert.empty()) return;
                for (size_t i = 0; i<m.vert.size(); ++i)
                {
                    if (pu.remap[i]<size_t(m.vn))
                    {
                        assert(!m.vert[i].IsD());
                        m.vert[pu.remap[i]].ImportData(m.vert[i]);
                        if (HasVFAdjacency(m))
                        {
                            if (m.vert[i].IsVFInitialized())
                            {
                                m.vert[pu.remap[i]].VFp() = m.vert[i].cVFp();
                                m.vert[pu.remap[i]].VFi() = m.vert[i].cVFi();
                            }
                            else m.vert[pu.remap[i]].VFClear();
                        }
                    }
                }

                // reorder the optional atttributes in m.vert_attr to reflect the changes
                ReorderAttribute(m.vert_attr, pu.remap, m);

                // setup the pointer updater
                pu.oldBase = &m.vert[0];
                pu.oldEnd = &m.vert.back() + 1;

                // resize
                m.vert.resize(m.vn);

                // setup the pointer updater
                pu.newBase = (m.vert.empty()) ? 0 : &m.vert[0];
                pu.newEnd = (m.vert.empty()) ? 0 : &m.vert.back() + 1;

                // resize the optional atttributes in m.vert_attr to reflect the changes
                ResizeAttribute(m.vert_attr, m.vn, m);

                // Loop on the face to update the pointers FV relation (vertex refs)
                for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
                    if (!(*fi).IsD())
                        for (int i = 0; i<fi->VN(); ++i)
                        {
                            size_t oldIndex = (*fi).V(i) - pu.oldBase;
                            assert(pu.oldBase <= (*fi).V(i) && oldIndex < pu.remap.size());
                            (*fi).V(i) = pu.newBase + pu.remap[oldIndex];
                        }
            }

        }; // end Allocator class


        template <class ComputeMeshType>
        class   UpdateNormal
        {
        public:
            typedef ComputeMeshType MeshType;
            typedef typename MeshType::VertexType     VertexType;
            typedef typename MeshType::CoordType     CoordType;
            typedef typename VertexType::NormalType     NormalType;
            typedef typename VertexType::ScalarType ScalarType;
            typedef typename MeshType::VertexPointer  VertexPointer;
            typedef typename MeshType::VertexIterator VertexIterator;
            typedef typename MeshType::FaceType       FaceType;
            typedef typename MeshType::FacePointer    FacePointer;
            typedef typename MeshType::FaceIterator   FaceIterator;

            /// \brief Set to zero all the PerVertex normals
            /**
            Set to zero all the PerVertex normals. Used by all the face averaging algorithms.
            by default it does not clear the normals of unreferenced vertices because they could be still useful
            */
            static void PerVertexClear(ComputeMeshType &m, bool ClearAllVertNormal = false)
            {
                RequirePerVertexNormal(m);
                if (ClearAllVertNormal)
                    UpdateFlags<ComputeMeshType>::VertexClearV(m);
                else
                {
                    UpdateFlags<ComputeMeshType>::VertexSetV(m);
                    for (FaceIterator f = m.face.begin(); f != m.face.end(); ++f)
                        if (!(*f).IsD())
                            for (int i = 0; i<3; ++i) (*f).V(i)->ClearV();
                }
                VertexIterator vi;
                for (vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                    if (!(*vi).IsD() && (*vi).IsRW() && (!(*vi).IsV()))
                        (*vi).N() = NormalType((ScalarType)0, (ScalarType)0, (ScalarType)0);
            }


            /// \brief Calculates the face normal
            ///
            /// Not normalized. Use PerFaceNormalized() or call NormalizePerVertex() if you need unit length per face normals.
            static void PerFace(ComputeMeshType &m)
            {
                RequirePerFaceNormal(m);
                for (FaceIterator f = m.face.begin(); f != m.face.end(); ++f)
                    if (!(*f).IsD())
                        f->N() = TriangleNormal(*f).Normalize();
            }


            static void PerVertexAngleWeighted(ComputeMeshType &m)
            {
                PerVertexClear(m);
                FaceIterator f;
                for (f = m.face.begin(); f != m.face.end(); ++f)
                    if (!(*f).IsD() && (*f).IsR())
                    {
                        NormalType t = TriangleNormal(*f).Normalize();
                        NormalType e0 = ((*f).V1(0)->cP() - (*f).V0(0)->cP()).Normalize();
                        NormalType e1 = ((*f).V1(1)->cP() - (*f).V0(1)->cP()).Normalize();
                        NormalType e2 = ((*f).V1(2)->cP() - (*f).V0(2)->cP()).Normalize();

                        (*f).V(0)->N() += t*AngleN(e0, -e2);
                        (*f).V(1)->N() += t*AngleN(-e0, e1);
                        (*f).V(2)->N() += t*AngleN(-e1, e2);
                    }
            } 

            /// \brief Equivalent to PerFace() and NormalizePerFace()
            static void PerFaceNormalized(ComputeMeshType &m)
            {
                PerFace(m);
                NormalizePerFace(m);
            }

            /// \brief Normalize the length of the face normals.
            static void NormalizePerFace(ComputeMeshType &m)
            {
                tri::RequirePerFaceNormal(m);
                for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
                    if (!(*fi).IsD())	
                        (*fi).N().Normalize();
            }
        }; // end class


        /** \brief Set the face incremental mark of the vertex to the one of the mesh.
        @param m the mesh containing the element
        @param f Vertex pointer */
        template <class MeshType> inline void Mark(MeshType & m, typename MeshType::FacePointer f) { f->IMark() = m.imark; }

        /** \brief Check if the face incremental mark matches the one of the mesh.
        @param m the mesh containing the element
        @param f Face pointer */
        template <class MeshType> inline bool IsMarked(MeshType & m, typename MeshType::ConstFacePointer f) { return f->cIMark() == m.imark; }

        /** \brief Unmark, in constant time, all the elements (face and vertices) of a mesh.
        @param m the mesh containing the element

        In practice this function just increment the internal counter that stores the value for which an element is considered marked;
        therefore all the mesh elements become immediately un-mmarked.
        */
        template <class MeshType> inline void UnMarkAll(MeshType & m)
        {
            ++m.imark;
        }


        /// Returns the normal to the plane passing through p0,p1,p2
        template<class TriangleType>
        static typename TriangleType::CoordType TriangleNormal(const TriangleType &t)
        {
            return ((t.cP(1) - t.cP(0)) ^ (t.cP(2) - t.cP(0)));
        }
    }	 // end namespace tri
}	 // end namespace

