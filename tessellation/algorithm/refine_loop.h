#pragma once

#include <cmath>
#include <algorithm/refine.h>

namespace tessellation
{
    namespace tri
    {

        /*
        Metodo di Loop dalla documentazione "Siggraph 2000 course on subdivision"

                d4------d3							d4------d3
             /	\		 / 	\						 /	\		 / 	\							u
            /		 \  /  	 \					/		e4--e3 	 \					 / \
         /	   	\/	 		\				 /	 / 	\/	\		\					/		\
        d5------d1------d2	->	d5--e5--d1--e2--d2			 l--M--r
         \	   	/\	 		/				 \	 \ 	/\	/		/				  \	  /
            \		 /  \ 	 /					\		e6--e7	 /					 \ /
             \	/		 \ 	/						 \	/		 \ 	/							d
                d6------d7							d6------d7

        *******************************************************

        */

        /*!
         * \brief Weight class for classical Loop's scheme.
         *
         * See Zorin, D. & Schr?eder, P.: Subdivision for modeling and animation. Proc. ACM SIGGRAPH [Courses], 2000
         */
        template<class SCALAR_TYPE>
        struct LoopWeight
        {
            inline SCALAR_TYPE beta(int k)
            {
                return (k > 3) ? (5.0 / 8.0 - std::pow((3.0 / 8.0 + std::cos(2.0 * M_PI / SCALAR_TYPE(k)) / 4.0), 2)) / SCALAR_TYPE(k) : 3.0 / 16.0;
            }

            inline SCALAR_TYPE incidentRegular(int)
            {
                return 3.0 / 8.0;
            }
            inline SCALAR_TYPE incidentIrregular(int)
            {
                return 3.0 / 8.0;
            }
            inline SCALAR_TYPE opposite(int)
            {
                return 1.0 / 8.0;
            }
        };

        // Centroid and LS3Projection classes may be pettre placed in an other file. (which one ?)

        /*!
         * \brief Allow to compute classical Loop subdivision surface with the same code than LS3.
         */
        template<class MESH_TYPE, class LSCALAR_TYPE = typename MESH_TYPE::ScalarType>
        struct Centroid
        {
            typedef typename MESH_TYPE::ScalarType Scalar;
            typedef typename MESH_TYPE::CoordType CoordType;
            typedef LSCALAR_TYPE LScalar;
            typedef tessellation::Point3<LScalar> LVector;

            LVector sumP;
            LScalar sumW;

            Centroid()
            {
                reset();
            }
            inline void reset()
            {
                sumP.SetZero();
                sumW = 0.;
            }
            inline void addVertex(const typename MESH_TYPE::VertexType &v, LScalar w)
            {
                LVector p(v.cP().X(), v.cP().Y(), v.cP().Z());
                sumP += p * w;
                sumW += w;
            }
            inline void project(std::pair<CoordType, CoordType> &nv) const
            {
                LVector position = sumP / sumW;
                nv.first = CoordType(position.X(), position.Y(), position.Z());
            }
        };

         template<class MESH_TYPE, class METHOD_TYPE = Centroid<MESH_TYPE>, class WEIGHT_TYPE = LoopWeight<typename MESH_TYPE::ScalarType> >
        struct OddPointLoopGeneric : public std::unary_function<face::Pos<typename MESH_TYPE::FaceType> , typename MESH_TYPE::VertexType>
        {
            typedef METHOD_TYPE Projection;
            typedef WEIGHT_TYPE Weight;
            typedef typename MESH_TYPE::template PerVertexAttributeHandle<int> ValenceAttr;
            typedef typename MESH_TYPE::CoordType CoordType;

            MESH_TYPE &m;
            Projection proj;
            Weight weight;
            ValenceAttr *valence;

            inline OddPointLoopGeneric(MESH_TYPE &_m, Projection proj = Projection(), Weight weight = Weight()) :
                m(_m), proj(proj), weight(weight), valence(0) {}

            void operator()(typename MESH_TYPE::VertexType &nv, face::Pos<typename MESH_TYPE::FaceType>  ep)
            {
                proj.reset();

                face::Pos<typename MESH_TYPE::FaceType> he(ep.f, ep.z, ep.f->V(ep.z));
                typename MESH_TYPE::VertexType *l, *r, *u, *d;
                l = he.v;
                he.FlipV();
                r = he.v;

                if( tri::HasPerVertexColor(m))
                    nv.C().lerp(ep.f->V(ep.z)->C(), ep.f->V1(ep.z)->C(), .5f);

                if (he.IsBorder())
                {
                    proj.addVertex(*l, 0.5);
                    proj.addVertex(*r, 0.5);
                    std::pair<CoordType, CoordType>pp;
                    proj.project(pp);
                    nv.P() = pp.first;
                    nv.N() = pp.second;
                }
                else
                {
                    he.FlipE();
                    he.FlipV();
                    u = he.v;
                    he.FlipV();
                    he.FlipE();
                    assert(he.v == r); // back to r
                    he.FlipF();
                    he.FlipE();
                    he.FlipV();
                    d = he.v;

                    if(valence && ((*valence)[l] == 6 || (*valence)[r] == 6))
                    {
                        int extra = ((*valence)[l] == 6) ? (*valence)[r] : (*valence)[l];
                        proj.addVertex(*l, ((*valence)[l] == 6) ? weight.incidentRegular(extra) : weight.incidentIrregular(extra));
                        proj.addVertex(*r, ((*valence)[r] == 6) ? weight.incidentRegular(extra) : weight.incidentIrregular(extra));
                        proj.addVertex(*u, weight.opposite(extra));
                        proj.addVertex(*d, weight.opposite(extra));
                    }
                    // unhandled case that append only at first subdivision step: use Loop's weights
                    else
                    {
                        proj.addVertex(*l, 3.0 / 8.0);
                        proj.addVertex(*r, 3.0 / 8.0);
                        proj.addVertex(*u, 1.0 / 8.0);
                        proj.addVertex(*d, 1.0 / 8.0);
                    }
                    std::pair<CoordType, CoordType>pp;
                    proj.project(pp);
                    nv.P() = pp.first;
                    nv.N() = pp.second;
                    //			proj.project(nv);
                }

            }

            template<class FL_TYPE>
            TexCoord2<FL_TYPE, 1> WedgeInterp(TexCoord2<FL_TYPE, 1> &t0, TexCoord2<FL_TYPE, 1> &t1)
            {
                TexCoord2<FL_TYPE, 1> tmp;
                tmp.n() = t0.n();
                tmp.t() = (t0.t() + t1.t()) / 2.0;
                return tmp;
            }

            inline void setValenceAttr(ValenceAttr *valence)
            {
                this->valence = valence;
            }
        };

        template<class MESH_TYPE, class METHOD_TYPE = Centroid<MESH_TYPE>, class WEIGHT_TYPE = LoopWeight<typename MESH_TYPE::ScalarType> >
        struct EvenPointLoopGeneric : public std::unary_function<face::Pos<typename MESH_TYPE::FaceType> , typename MESH_TYPE::VertexType>
        {
            typedef METHOD_TYPE Projection;
            typedef WEIGHT_TYPE Weight;
            typedef typename MESH_TYPE::template PerVertexAttributeHandle<int> ValenceAttr;
            typedef typename MESH_TYPE::CoordType CoordType;

            Projection proj;
            Weight weight;
            ValenceAttr *valence;

            inline EvenPointLoopGeneric(Projection proj = Projection(), Weight weight = Weight()) :
                proj(proj), weight(weight), valence(0) {}

            void operator()(std::pair<CoordType, CoordType> &nv, face::Pos<typename MESH_TYPE::FaceType>  ep)
            {
                proj.reset();

                face::Pos<typename MESH_TYPE::FaceType> he(ep.f, ep.z, ep.f->V(ep.z));
                typename MESH_TYPE::VertexType *r, *l,  *curr;
                curr = he.v;
                face::Pos<typename MESH_TYPE::FaceType> heStart = he;

                // compute valence of this vertex or find a border
                int k = 0;
                do
                {
                    he.NextE();
                    k++;
                }
                while(!he.IsBorder() && he != heStart);

                if (he.IsBorder())  	//	Border rule
                {
                    // consider valence of borders as if they are half+1 of an inner vertex. not perfect, but better than nothing.
                    if(valence)
                    {
                        k = 0;
                        do
                        {
                            he.NextE();
                            k++;
                        }
                        while(!he.IsBorder());
                        (*valence)[he.V()] = std::max(2 * (k - 1), 3);
                        //				(*valence)[he.V()] = 6;
                    }

                    he.FlipV();
                    r = he.v;
                    he.FlipV();
                    he.NextB();
                    l = he.v;

                    proj.addVertex(*curr, 3.0 / 4.0);
                    proj.addVertex(*l, 1.0 / 8.0);
                    proj.addVertex(*r, 1.0 / 8.0);
                    //			std::pair<Point3f,Point3f>pp;
                    proj.project(nv);
                    //			nv.P()=pp.first;
                    //			nv.N()=pp.second;
                    //			proj.project(nv);
                }
                else  	//	Inner rule
                {
                    //			assert(!he.v->IsB()); border flag no longer updated (useless)
                    if(valence)
                        (*valence)[he.V()] = k;

                    typename MESH_TYPE::ScalarType beta = weight.beta(k);

                    proj.addVertex(*curr, 1.0 - (typename MESH_TYPE::ScalarType)(k) * beta);
                    do
                    {
                        proj.addVertex(*he.VFlip(), beta);
                        he.NextE();
                    }
                    while(he != heStart);

                    proj.project(nv);
                }
            } // end of operator()


            inline void setValenceAttr(ValenceAttr *valence)
            {
                this->valence = valence;
            }
        };

        template<class MESH_TYPE>
        struct OddPointLoop : OddPointLoopGeneric<MESH_TYPE, Centroid<MESH_TYPE> >
        {
            OddPointLoop(MESH_TYPE &_m) : OddPointLoopGeneric<MESH_TYPE, Centroid<MESH_TYPE> >(_m) {}
        };

        template<class MESH_TYPE>
        struct EvenPointLoop : EvenPointLoopGeneric<MESH_TYPE, Centroid<MESH_TYPE> >
        {
        };

        template<class MESH_TYPE, class ODD_VERT, class EVEN_VERT>
        bool RefineOddEven(MESH_TYPE &m, ODD_VERT odd, EVEN_VERT even, float length,
                           bool RefineSelected = false)
        {
            EdgeLen <MESH_TYPE, typename MESH_TYPE::ScalarType> ep(length);
            return RefineOddEvenE(m, odd, even, ep, RefineSelected);
        }

        /*!
         * \brief Perform diadic subdivision using given rules for odd and even vertices.
         */
        template<class MESH_TYPE, class ODD_VERT, class EVEN_VERT, class PREDICATE>
        bool RefineOddEvenE(MESH_TYPE &m, ODD_VERT odd, EVEN_VERT even, PREDICATE edgePred,
                            bool RefineSelected = false)
        {
            typedef typename MESH_TYPE::template PerVertexAttributeHandle<int> ValenceAttr;

            // momentaneamente le callback sono identiche, almeno cbOdd deve essere passata
            //cbEven = cbOdd;

            // to mark visited vertices
            int evenFlag = MESH_TYPE::VertexType::NewBitFlag();
            for (int i = 0; i < m.vn ; i++ )
            {
                m.vert[i].ClearUserBit(evenFlag);
            }

            ValenceAttr valence = tessellation::tri::Allocator<MESH_TYPE>:: template AddPerVertexAttribute<int>(m);
            odd.setValenceAttr(&valence);
            even.setValenceAttr(&valence);

            // store updated vertices
            std::vector<bool> updatedList(m.vn, false);
            //std::vector<typename MESH_TYPE::VertexType> newEven(m.vn);
            std::vector<std::pair<typename MESH_TYPE::CoordType, typename MESH_TYPE::CoordType> > newEven(m.vn);

            typename MESH_TYPE::VertexIterator vi;
            typename MESH_TYPE::FaceIterator fi;
            for (fi = m.face.begin(); fi != m.face.end(); fi++) if(!(*fi).IsD() && (!RefineSelected || (*fi).IsS()))  //itero facce
                {
                    for (int i = 0; i < 3; i++)   //itero vert
                    {
                        if ( !(*fi).V(i)->IsUserBit(evenFlag) && ! (*fi).V(i)->IsD() )
                        {
                            (*fi).V(i)->SetUserBit(evenFlag);
                            //	use face selection, not vertex selection, to be coherent with RefineE
                            //if (RefineSelected && !(*fi).V(i)->IsS() )
                            //	break;
                            face::Pos<typename MESH_TYPE::FaceType>aux (&(*fi), i);
                            if( tri::HasPerVertexColor(m) )
                            {
                                (*fi).V(i)->C().lerp((*fi).V0(i)->C() , (*fi).V1(i)->C(), 0.5f);
                            }

                            //if (cbEven) {
                            //    (*cbEven)(int(100.0f * (float)j / (float)m.fn),"Refining");
                            //    j++;
                            //}
                            size_t index = tri::Index(m, (*fi).V(i));
                            updatedList[index] = true;
                            even(newEven[index], aux);
                        }
                    }
                }

            MESH_TYPE::VertexType::DeleteBitFlag(evenFlag);

            // Now apply the stored normal and position to the initial vertex set (note that newEven is << m.vert)
            RefineE< MESH_TYPE, ODD_VERT > (m, odd, edgePred, RefineSelected);
            for(size_t i = 0; i < newEven.size(); ++i)
            {
                if(updatedList[i])
                {
                    m.vert[i].P() = newEven[i].first;
                    m.vert[i].N() = newEven[i].second;
                }
            }

            odd.setValenceAttr(0);
            even.setValenceAttr(0);

            tessellation::tri::Allocator<MESH_TYPE>::DeletePerVertexAttribute(m, valence);

            return true;
        }

    } // namespace tri
} // namespace tessellation

