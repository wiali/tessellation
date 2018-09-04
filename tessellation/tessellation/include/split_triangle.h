#pragma once

#include <cmath>

namespace tessellation
{
    namespace tri
    {
        struct Split
        {
            int TriNum;
            int TV[4][3];
            int swap[2][2];
            int TE[4][3];
        };

        const Split SplitTab[8] =
        {
            { 1, { { 0, 1, 2 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }, { { 0, 0 }, { 0, 0 } }, { { 0, 1, 2 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { 2, { { 0, 3, 2 }, { 3, 1, 2 }, { 0, 0, 0 }, { 0, 0, 0 } }, { { 0, 0 }, { 0, 0 } }, { { 0, 3, 2 }, { 0, 1, 3 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { 2, { { 0, 1, 4 }, { 0, 4, 2 }, { 0, 0, 0 }, { 0, 0, 0 } }, { { 0, 0 }, { 0, 0 } }, { { 0, 1, 3 }, { 3, 1, 2 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { 3, { { 3, 1, 4 }, { 0, 3, 2 }, { 4, 2, 3 }, { 0, 0, 0 } }, { { 0, 4 }, { 3, 2 } }, { { 0, 1, 3 }, { 0, 3, 2 }, { 1, 3, 3 }, { 0, 0, 0 } } },
            { 2, { { 0, 1, 5 }, { 5, 1, 2 }, { 0, 0, 0 }, { 0, 0, 0 } }, { { 0, 0 }, { 0, 0 } }, { { 0, 3, 2 }, { 3, 1, 2 }, { 0, 0, 0 }, { 0, 0, 0 } } },
            { 3, { { 0, 3, 5 }, { 3, 1, 5 }, { 2, 5, 1 }, { 0, 0, 0 } }, { { 3, 2 }, { 5, 1 } }, { { 0, 3, 2 }, { 0, 3, 3 }, { 2, 3, 1 }, { 0, 0, 0 } } },
            { 3, { { 2, 5, 4 }, { 0, 1, 5 }, { 4, 5, 1 }, { 0, 0, 0 } }, { { 0, 4 }, { 5, 1 } }, { { 2, 3, 1 }, { 0, 3, 2 }, { 3, 3, 1 }, { 0, 0, 0 } } },
            { 4, { { 3, 4, 5 }, { 0, 3, 5 }, { 3, 1, 4 }, { 5, 4, 2 } }, { { 0, 0 }, { 0, 0 } }, { { 3, 3, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 3, 1, 2 } } },
        };

        template <class MESH_TYPE, class FLT>
        class EdgeLen
        {
            FLT squaredThr;
        public:
            EdgeLen() {};
            EdgeLen(FLT threshold)
            {
                setThr(threshold);
            }
            void setThr(FLT threshold)
            {
                squaredThr = threshold * threshold;
            }
            bool operator()(face::Pos<typename MESH_TYPE::FaceType> ep) const
            {
                return SquaredDistance(ep.V()->P(), ep.VFlip()->P()) > squaredThr;
            }
        };

        template <class VertexPointer>
        class RefinedFaceData
        {
        public:
            RefinedFaceData()
            {
                ep[0] = 0;
                ep[1] = 0;
                ep[2] = 0;
                vp[0] = 0;
                vp[1] = 0;
                vp[2] = 0;
            }
            bool ep[3];
            VertexPointer vp[3];
        };

        template<class MESH_TYPE, class MIDPOINT, class EDGEPRED>
        bool RefineE(MESH_TYPE& m, MIDPOINT& mid, EDGEPRED& ep, bool RefineSelected = false)
        {
            typedef typename MESH_TYPE::VertexIterator VertexIterator;
            typedef typename MESH_TYPE::FaceIterator FaceIterator;
            typedef typename MESH_TYPE::VertexPointer VertexPointer;
            typedef typename MESH_TYPE::FacePointer FacePointer;
            typedef typename MESH_TYPE::FaceType FaceType;
            typedef typename MESH_TYPE::FaceType::TexCoordType TexCoordType;
            assert(tri::HasFFAdjacency(m));
            tri::UpdateFlags<MESH_TYPE>::FaceBorderFromFF(m);
            typedef face::Pos<FaceType> PosType;
            int j, NewVertNum = 0, NewFaceNum = 0;
            typedef RefinedFaceData<VertexPointer> RFD;
            typedef typename MESH_TYPE :: template PerFaceAttributeHandle<RFD> HandleType;
            HandleType RD = tri::Allocator<MESH_TYPE>:: template AddPerFaceAttribute<RFD>(m, std::string("RefineData"));
            FaceIterator fi;
            for (fi = m.face.begin(), j = 0; fi != m.face.end(); ++fi) if (!(*fi).IsD())
                {
                    if (RefineSelected && !(*fi).IsS())
                    {
                        continue;
                    }
                    for (j = 0; j < 3; j++)
                    {
                        if (RD[fi].ep[j])
                        {
                            continue;
                        }
                        PosType edgeCur(&*fi, j);
                        if (RefineSelected && !edgeCur.FFlip()->IsS())
                        {
                            continue;
                        }
                        if (!ep(edgeCur))
                        {
                            continue;
                        }
                        RD[edgeCur.F()].ep[edgeCur.E()] = true;
                        ++NewFaceNum;
                        ++NewVertNum;
                        if (!edgeCur.IsBorder())
                        {
                            edgeCur.FlipF();
                            edgeCur.F()->SetV();
                            RD[edgeCur.F()].ep[edgeCur.E()] = true;
                            ++NewFaceNum;
                        }
                    }
                }
            if (NewVertNum == 0)
            {
                tri::Allocator<MESH_TYPE>::template DeletePerFaceAttribute<RefinedFaceData<VertexPointer>>(m, RD);
                return false;
            }
            VertexIterator lastv = tri::Allocator<MESH_TYPE>::AddVertices(m, NewVertNum);
            for (fi = m.face.begin(); fi != m.face.end(); ++fi) if (!(*fi).IsD())
            {
                for (j = 0; j < 3; j++)
                {
                    if (RefineSelected && !(*fi).IsS())
                    {
                        continue;
                    }
                    for (j = 0; j < 3; j++)
                    {
                        PosType edgeCur(&*fi, j);
                        if (RefineSelected && !edgeCur.FFlip()->IsS())
                        {
                            continue;
                        }
                        if (RD[edgeCur.F()].ep[edgeCur.E()] && RD[edgeCur.F()].vp[edgeCur.E()] == 0)
                        {
                            RD[edgeCur.F()].vp[edgeCur.E()] = &*lastv;
                            mid(*lastv, edgeCur);
                            if (!edgeCur.IsBorder())
                            {
                                edgeCur.FlipF();
                                assert(RD[edgeCur.F()].ep[edgeCur.E()]);
                                RD[edgeCur.F()].vp[edgeCur.E()] = &*lastv;
                            }
                            ++lastv;
                        }
                    }
                }
            }
            assert(lastv == m.vert.end());
            FaceIterator lastf = tri::Allocator<MESH_TYPE>::AddFaces(m, NewFaceNum);
            FaceIterator oldendf = lastf;
            VertexPointer vv[6];
            FacePointer nf[4];
            TexCoordType wtt[6];
            int fca = 0, fcn = 0;
            for (fi = m.face.begin(); fi != oldendf; ++fi) if (!(*fi).IsD())
            {
                fcn++;
                vv[0] = (*fi).V(0);
                vv[1] = (*fi).V(1);
                vv[2] = (*fi).V(2);
                vv[3] = RD[fi].vp[0];
                vv[4] = RD[fi].vp[1];
                vv[5] = RD[fi].vp[2];
                int ind = ((&*vv[3]) ? 1 : 0) + ((&*vv[4]) ? 2 : 0) + ((&*vv[5]) ? 4 : 0);
                nf[0] = &*fi;
                int i;
                for (i = 1; i < SplitTab[ind].TriNum; ++i)
                {
                    nf[i] = &*lastf;
                    ++lastf;
                    fca++;
                    if (RefineSelected || (*fi).IsS())
                    {
                        (*nf[i]).SetS();
                    }
                    nf[i]->ImportData(*fi);
                }
                int orgflag = (*fi).Flags();
                for (i = 0; i < SplitTab[ind].TriNum; ++i)
                    for (j = 0; j < 3; ++j)
                    {
                        (*nf[i]).V(j) = &*vv[SplitTab[ind].TV[i][j]];
                        assert((*nf[i]).V(j) != 0);
                        if (SplitTab[ind].TE[i][j] != 3)
                        {
                            if (orgflag & (MESH_TYPE::FaceType::BORDER0 << (SplitTab[ind].TE[i][j])))
                            {
                                (*nf[i]).SetB(j);
                            }
                            else
                            {
                                (*nf[i]).ClearB(j);
                            }
                        }
                        else
                        {
                            (*nf[i]).ClearB(j);
                        }
                    }
                if (SplitTab[ind].TriNum == 3 &&
                        SquaredDistance(vv[SplitTab[ind].swap[0][0]]->P(), vv[SplitTab[ind].swap[0][1]]->P()) <
                        SquaredDistance(vv[SplitTab[ind].swap[1][0]]->P(), vv[SplitTab[ind].swap[1][1]]->P()))
                {
                    (*nf[2]).V(1) = (*nf[1]).V(0);
                    (*nf[1]).V(1) = (*nf[2]).V(0);
                    if ((*nf[1]).IsB(0))
                    {
                        (*nf[2]).SetB(1);
                    }
                    else
                    {
                        (*nf[2]).ClearB(1);
                    }
                    if ((*nf[2]).IsB(0))
                    {
                        (*nf[1]).SetB(1);
                    }
                    else
                    {
                        (*nf[1]).ClearB(1);
                    }
                    (*nf[1]).ClearB(0);
                    (*nf[2]).ClearB(0);
                }
            }
            assert(lastf == m.face.end());
            assert(!m.vert.empty());
            for (fi = m.face.begin(); fi != m.face.end(); ++fi) if (!(*fi).IsD())
            {
                assert((*fi).V(0) >= &*m.vert.begin() && (*fi).V(0) <= &m.vert.back());
                assert((*fi).V(1) >= &*m.vert.begin() && (*fi).V(1) <= &m.vert.back());
                assert((*fi).V(2) >= &*m.vert.begin() && (*fi).V(2) <= &m.vert.back());
            }
            tri::UpdateTopology<MESH_TYPE>::FaceFace(m);
            tri::Allocator<MESH_TYPE> :: template DeletePerFaceAttribute<RefinedFaceData<VertexPointer> >(m, RD);
            return true;
        }

        template<class MESH_TYPE, class MIDPOINT>
        bool Refine(MESH_TYPE& m, MIDPOINT mid, typename MESH_TYPE::ScalarType thr = 0, bool RefineSelected = false)
        {
            EdgeLen <MESH_TYPE, typename MESH_TYPE::ScalarType> ep(thr);
            return RefineE(m, mid, ep, RefineSelected);
        }

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
            inline void addVertex(const typename MESH_TYPE::VertexType& v, LScalar w)
            {
                LVector p(v.cP().X(), v.cP().Y(), v.cP().Z());
                sumP += p * w;
                sumW += w;
            }
            inline void project(std::pair<CoordType, CoordType>& nv) const
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

            MESH_TYPE& m;
            Projection proj;
            Weight weight;
            ValenceAttr* valence;

            inline OddPointLoopGeneric(MESH_TYPE& _m, Projection proj = Projection(), Weight weight = Weight()) :
                m(_m), proj(proj), weight(weight), valence(0) {}

            void operator()(typename MESH_TYPE::VertexType& nv, face::Pos<typename MESH_TYPE::FaceType>  ep)
            {
                proj.reset();
                face::Pos<typename MESH_TYPE::FaceType> he(ep.f, ep.z, ep.f->V(ep.z));
                typename MESH_TYPE::VertexType* l, *r, *u, *d;
                l = he.v;
                he.FlipV();
                r = he.v;
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
                    assert(he.v == r);
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
                }
            }

            inline void setValenceAttr(ValenceAttr* valence)
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
            ValenceAttr* valence;

            inline EvenPointLoopGeneric(Projection proj = Projection(), Weight weight = Weight()) :
                proj(proj), weight(weight), valence(0) {}

            void operator()(std::pair<CoordType, CoordType>& nv, face::Pos<typename MESH_TYPE::FaceType>  ep)
            {
                proj.reset();
                face::Pos<typename MESH_TYPE::FaceType> he(ep.f, ep.z, ep.f->V(ep.z));
                typename MESH_TYPE::VertexType* r, *l,  *curr;
                curr = he.v;
                face::Pos<typename MESH_TYPE::FaceType> heStart = he;
                int k = 0;
                do
                {
                    he.NextE();
                    k++;
                }
                while(!he.IsBorder() && he != heStart);
                if (he.IsBorder())
                {
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
                    }
                    he.FlipV();
                    r = he.v;
                    he.FlipV();
                    he.NextB();
                    l = he.v;
                    proj.addVertex(*curr, 3.0 / 4.0);
                    proj.addVertex(*l, 1.0 / 8.0);
                    proj.addVertex(*r, 1.0 / 8.0);
                    proj.project(nv);
                }
                else
                {
                    if(valence)
                    {
                        (*valence)[he.V()] = k;
                    }
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
            }

            inline void setValenceAttr(ValenceAttr* valence)
            {
                this->valence = valence;
            }
        };

        template<class MESH_TYPE>
        struct OddPointLoop : OddPointLoopGeneric<MESH_TYPE, Centroid<MESH_TYPE> >
        {
            OddPointLoop(MESH_TYPE& _m) : OddPointLoopGeneric<MESH_TYPE, Centroid<MESH_TYPE> >(_m) {}
        };

        template<class MESH_TYPE>
        struct EvenPointLoop : EvenPointLoopGeneric<MESH_TYPE, Centroid<MESH_TYPE> >
        {
        };

        template<class MESH_TYPE, class ODD_VERT, class EVEN_VERT>
        bool RefineOddEven(MESH_TYPE& m, ODD_VERT odd, EVEN_VERT even, float length,
                           bool RefineSelected = false)
        {
            EdgeLen <MESH_TYPE, typename MESH_TYPE::ScalarType> ep(length);
            return RefineOddEvenE(m, odd, even, ep, RefineSelected);
        }

        template<class MESH_TYPE, class ODD_VERT, class EVEN_VERT, class PREDICATE>
        bool RefineOddEvenE(MESH_TYPE& m, ODD_VERT odd, EVEN_VERT even, PREDICATE edgePred,
                            bool RefineSelected = false)
        {
            typedef typename MESH_TYPE::template PerVertexAttributeHandle<int> ValenceAttr;
            int evenFlag = MESH_TYPE::VertexType::NewBitFlag();
            for (int i = 0; i < m.vn ; i++ )
            {
                m.vert[i].ClearUserBit(evenFlag);
            }
            ValenceAttr valence = tessellation::tri::Allocator<MESH_TYPE>:: template AddPerVertexAttribute<int>(m);
            odd.setValenceAttr(&valence);
            even.setValenceAttr(&valence);
            std::vector<bool> updatedList(m.vn, false);
            std::vector<std::pair<typename MESH_TYPE::CoordType, typename MESH_TYPE::CoordType> > newEven(m.vn);
            typename MESH_TYPE::VertexIterator vi;
            typename MESH_TYPE::FaceIterator fi;
            for (fi = m.face.begin(); fi != m.face.end(); fi++) if(!(*fi).IsD() && (!RefineSelected || (*fi).IsS()))
                {
                    for (int i = 0; i < 3; i++)
                    {
                        if ( !(*fi).V(i)->IsUserBit(evenFlag) && ! (*fi).V(i)->IsD() )
                        {
                            (*fi).V(i)->SetUserBit(evenFlag);
                            face::Pos<typename MESH_TYPE::FaceType>aux (&(*fi), i);
                            size_t index = tri::Index(m, (*fi).V(i));
                            updatedList[index] = true;
                            even(newEven[index], aux);
                        }
                    }
                }
            MESH_TYPE::VertexType::DeleteBitFlag(evenFlag);
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
    }
}

