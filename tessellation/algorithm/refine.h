#pragma once

#include <functional>
#include <map>

namespace tessellation
{
    namespace tri
    {

        /* A very short intro about the generic refinement framework,
            the main fuction is the

         template<class MESH_TYPE,class MIDPOINT, class EDGEPRED>
         bool RefineE(MESH_TYPE &m, MIDPOINT mid, EDGEPRED ep,bool RefineSelected=false, CallBackPos *cb = 0)

         You have to provide two functor objects to this, one for deciding what edge has to be spltted and one to decide position and new values for the attributes of the new point.

         for example the minimal EDGEPRED is

         template <class MESH_TYPE, class FLT> class EdgeLen
         {
           public:
             FLT thr2;
             bool operator()(face::Pos<typename MESH_TYPE::FaceType> ep) const
             {
                    return SquaredDistance(ep.f->V(ep.z)->P(), ep.f->V1(ep.z)->P())>thr2;
             }
         };

         With a bit of patience you can customize to make also slicing operation.

        */


        /* The table which encodes how to subdivide a triangle depending
           on the splitted edges is organized as such:

            TriNum (the first number):    encodes the number of triangles
            TV (the following 4 triples): encodes the resulting triangles where
                  0, 1, 2 are the original vertices of the triangles and 3, 4, 5
                  (mp01, mp12, mp20) are the midpoints of the three edges.

           In the case two edges are splitted the triangle has 2 possible splittings:
        we need to choose a diagonal of the resulting trapezoid.
        'swap' encodes the two diagonals to test: if diag1 < diag2 we swap the diagonal
        like this (140, 504 -> 150, 514) (the second vertex of each triangles is replaced
         by the first vertex of the other one).
                    2
                   / \
                  5---4
                 /     \
                0-------1

        */

        class Split
        {
        public:
            int TriNum;			// number of triangles
            int TV[4][3];   // The triangles coded as the following convention
            //     0..2 vertici originali del triangolo
            //     3..5 mp01, mp12, mp20 midpoints of the three edges
            int swap[2][2]; // the two diagonals to test for swapping
            int TE[4][3];   // the edge-edge correspondence between refined triangles and the old one
            //      (3) means the edge of the new triangle is internal;
        };

        const Split SplitTab[8] =
        {
            /* m20 m12 m01 */
            /*  0   0   0 */	{1, {{0, 1, 2}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, 0}, {0, 0}},  {{0, 1, 2}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}} },
            /*  0   0   1 */	{2, {{0, 3, 2}, {3, 1, 2}, {0, 0, 0}, {0, 0, 0}}, {{0, 0}, {0, 0}},  {{0, 3, 2}, {0, 1, 3}, {0, 0, 0}, {0, 0, 0}} },
            /*  0   1   0 */	{2, {{0, 1, 4}, {0, 4, 2}, {0, 0, 0}, {0, 0, 0}}, {{0, 0}, {0, 0}},  {{0, 1, 3}, {3, 1, 2}, {0, 0, 0}, {0, 0, 0}} },
            /*  0   1   1 */	{3, {{3, 1, 4}, {0, 3, 2}, {4, 2, 3}, {0, 0, 0}}, {{0, 4}, {3, 2}},  {{0, 1, 3}, {0, 3, 2}, {1, 3, 3}, {0, 0, 0}} },
            /*  1   0   0 */	{2, {{0, 1, 5}, {5, 1, 2}, {0, 0, 0}, {0, 0, 0}}, {{0, 0}, {0, 0}},  {{0, 3, 2}, {3, 1, 2}, {0, 0, 0}, {0, 0, 0}} },
            /*  1   0   1 */	{3, {{0, 3, 5}, {3, 1, 5}, {2, 5, 1}, {0, 0, 0}}, {{3, 2}, {5, 1}},  {{0, 3, 2}, {0, 3, 3}, {2, 3, 1}, {0, 0, 0}} },
            /*  1   1   0 */	{3, {{2, 5, 4}, {0, 1, 5}, {4, 5, 1}, {0, 0, 0}}, {{0, 4}, {5, 1}},  {{2, 3, 1}, {0, 3, 2}, {3, 3, 1}, {0, 0, 0}} },
            /*  1   1   1 */	//{4, {{0,3,5},{3,1,4},{5,4,2},{3,4,5}}, {{0,0},{0,0}},  {{0,3,2},{0,1,3},{3,1,2},{3,3,3}} },
            /*  1   1   1 */	{4, {{3, 4, 5}, {0, 3, 5}, {3, 1, 4}, {5, 4, 2}}, {{0, 0}, {0, 0}},  {{3, 3, 3}, {0, 3, 2}, {0, 1, 3}, {3, 1, 2}} },
        };

        // Basic subdivision class
        // This class must provide methods for finding the position of the newly created vertices
        // In this implemenation we simply put the new vertex in the MidPoint position.
        // Color and TexCoords are interpolated accordingly.
        // This subdivision class allow also the correct interpolation of userdefined data by
        // providing, in the constructor, an interpolator functor that will be called for each new vertex to be created.

        

        // Basic Predicate that tells if a given edge must be splitted.
        // the constructure requires the threshold.
        // VERY IMPORTANT REQUIREMENT: this function must be symmetric
        // e.g. it must return the same value if the Pos is VFlipped.
        // If this function is not symmetric the Refine can crash.

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

        /*********************************************************/
        /*********************************************************

        Given a mesh the following function refines it according to two functor objects:

        - a predicate that tells if a given edge must be splitted

        - a functor that gives you the new poistion of the created vertices (starting from an edge)

        If RefineSelected is true only selected faces are taken into account for being splitted.

        Requirement: FF Adjacency and Manifoldness

        **********************************************************/
        /*********************************************************/
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
        bool RefineE(MESH_TYPE &m, MIDPOINT &mid, EDGEPRED &ep, bool RefineSelected = false)
        {
            // common typenames
            typedef typename MESH_TYPE::VertexIterator VertexIterator;
            typedef typename MESH_TYPE::FaceIterator FaceIterator;
            typedef typename MESH_TYPE::VertexPointer VertexPointer;
            typedef typename MESH_TYPE::FacePointer FacePointer;
            typedef typename MESH_TYPE::FaceType FaceType;
            typedef typename MESH_TYPE::FaceType::TexCoordType TexCoordType;
            assert(tri::HasFFAdjacency(m));
            tri::UpdateFlags<MESH_TYPE>::FaceBorderFromFF(m);
            typedef face::Pos<FaceType>  PosType;

            int j, NewVertNum = 0, NewFaceNum = 0;

            typedef RefinedFaceData<VertexPointer> RFD;
            typedef typename MESH_TYPE :: template PerFaceAttributeHandle<RFD> HandleType;
            HandleType RD  = tri::Allocator<MESH_TYPE>:: template AddPerFaceAttribute<RFD> (m, std::string("RefineData"));

            // First Loop: We analyze the mesh to compute the number of the new faces and new vertices
            FaceIterator fi;
            for(fi = m.face.begin(), j = 0; fi != m.face.end(); ++fi) if(!(*fi).IsD())
                {
                    // skip unselected faces if necessary
                    if(RefineSelected && !(*fi).IsS()) continue;

                    for(j = 0; j < 3; j++)
                    {
                        if(RD[fi].ep[j]) continue;

                        PosType edgeCur(&*fi, j);
                        if(RefineSelected && ! edgeCur.FFlip()->IsS()) continue;
                        if(!ep(edgeCur)) continue;

                        RD[edgeCur.F()].ep[edgeCur.E()] = true;
                        ++NewFaceNum;
                        ++NewVertNum;
                        //assert(edgeCur.IsManifold());
                        if(!edgeCur.IsBorder())
                        {
                            edgeCur.FlipF();
                            edgeCur.F()->SetV();
                            RD[edgeCur.F()].ep[edgeCur.E()] = true;
                            ++NewFaceNum;
                        }
                    }

                } // end face loop

            if(NewVertNum == 0 )
            {
                tri::Allocator<MESH_TYPE> :: template DeletePerFaceAttribute<RefinedFaceData<VertexPointer> >  (m, RD);
                return false;
            }
            VertexIterator lastv = tri::Allocator<MESH_TYPE>::AddVertices(m, NewVertNum);

            // Secondo loop: We initialize a edge->vertex map

            for(fi = m.face.begin(); fi != m.face.end(); ++fi) if(!(*fi).IsD())
                {
                    for(j = 0; j < 3; j++)
                    {
                        // skip unselected faces if necessary
                        if(RefineSelected && !(*fi).IsS()) continue;
                        for(j = 0; j < 3; j++)
                        {
                            PosType edgeCur(&*fi, j);
                            if(RefineSelected && ! edgeCur.FFlip()->IsS()) continue;

                            if( RD[edgeCur.F()].ep[edgeCur.E()] &&  RD[edgeCur.F()].vp[edgeCur.E()] == 0 )
                            {
                                RD[edgeCur.F()].vp[edgeCur.E()] = &*lastv;
                                mid(*lastv, edgeCur);
                                if(!edgeCur.IsBorder())
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

            assert(lastv == m.vert.end()); // critical assert: we MUST have used all the vertex that we forecasted we need

            FaceIterator lastf = tri::Allocator<MESH_TYPE>::AddFaces(m, NewFaceNum);
            FaceIterator oldendf = lastf;

            /*
             *               v0
             *
             *
             *               f0
             *
             *       mp01     f3     mp02
             *
             *
             *       f1               f2
             *
             *v1            mp12                v2
             *
            */

            VertexPointer vv[6];	// The six vertices that arise in the single triangle splitting
            //     0..2 Original triangle vertices
            //     3..5 mp01, mp12, mp20 midpoints of the three edges
            FacePointer nf[4];   // The (up to) four faces that are created.

            TexCoordType wtt[6];  // per ogni faccia sono al piu' tre i nuovi valori
            // di texture per wedge (uno per ogni edge)

            int fca = 0, fcn = 0;
            for(fi = m.face.begin(); fi != oldendf; ++fi) if(!(*fi).IsD())
                {
                    //if(cb && (++step%PercStep)==0)(*cb)(step/PercStep,"Refining...");
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
                    for(i = 1; i < SplitTab[ind].TriNum; ++i)
                    {
                        nf[i] = &*lastf;
                        ++lastf;
                        fca++;
                        if(RefineSelected || (*fi).IsS()) (*nf[i]).SetS();
                        nf[i]->ImportData(*fi);
                        //		if(tri::HasPerFaceColor(m))
                        //  		nf[i]->C()=(*fi).cC();
                    }


                    if(tri::HasPerWedgeTexCoord(m))
                        for(i = 0; i < 3; ++i)
                        {
                            wtt[i] = (*fi).WT(i);
                            wtt[3 + i] = mid.WedgeInterp((*fi).WT(i), (*fi).WT((i + 1) % 3));
                        }

                    int orgflag =	(*fi).Flags();
                    for(i = 0; i < SplitTab[ind].TriNum; ++i)
                        for(j = 0; j < 3; ++j)
                        {
                            (*nf[i]).V(j) = &*vv[SplitTab[ind].TV[i][j]];

                            if(tri::HasPerWedgeTexCoord(m)) //analogo ai vertici...
                                (*nf[i]).WT(j) = wtt[SplitTab[ind].TV[i][j]];

                            assert((*nf[i]).V(j) != 0);
                            if(SplitTab[ind].TE[i][j] != 3)
                            {
                                if(orgflag & (MESH_TYPE::FaceType::BORDER0 << (SplitTab[ind].TE[i][j])))
                                    (*nf[i]).SetB(j);
                                else
                                    (*nf[i]).ClearB(j);
                            }
                            else (*nf[i]).ClearB(j);
                        }

                    if(SplitTab[ind].TriNum == 3 &&
                            SquaredDistance(vv[SplitTab[ind].swap[0][0]]->P(), vv[SplitTab[ind].swap[0][1]]->P()) <
                            SquaredDistance(vv[SplitTab[ind].swap[1][0]]->P(), vv[SplitTab[ind].swap[1][1]]->P()) )
                    {
                        // swap the last two triangles
                        (*nf[2]).V(1) = (*nf[1]).V(0);
                        (*nf[1]).V(1) = (*nf[2]).V(0);
                        if(tri::HasPerWedgeTexCoord(m))  //swap also textures coordinates
                        {
                            (*nf[2]).WT(1) = (*nf[1]).WT(0);
                            (*nf[1]).WT(1) = (*nf[2]).WT(0);
                        }

                        if((*nf[1]).IsB(0)) (*nf[2]).SetB(1);
                        else (*nf[2]).ClearB(1);
                        if((*nf[2]).IsB(0)) (*nf[1]).SetB(1);
                        else (*nf[1]).ClearB(1);
                        (*nf[1]).ClearB(0);
                        (*nf[2]).ClearB(0);
                    }
                }

            assert(lastf == m.face.end());	 // critical assert: we MUST have used all the faces that we forecasted we need and that we previously allocated.
            assert(!m.vert.empty());
            for(fi = m.face.begin(); fi != m.face.end(); ++fi) if(!(*fi).IsD())
                {
                    assert((*fi).V(0) >= &*m.vert.begin() && (*fi).V(0) <= &m.vert.back() );
                    assert((*fi).V(1) >= &*m.vert.begin() && (*fi).V(1) <= &m.vert.back() );
                    assert((*fi).V(2) >= &*m.vert.begin() && (*fi).V(2) <= &m.vert.back() );
                }
            tri::UpdateTopology<MESH_TYPE>::FaceFace(m);

            tri::Allocator<MESH_TYPE> :: template DeletePerFaceAttribute<RefinedFaceData<VertexPointer> >  (m, RD);

            return true;
        }

        /*************************************************************************/
        // simple wrapper of the base refine for lazy coder that do not need a edge predicate

        template<class MESH_TYPE, class MIDPOINT>
        bool Refine(MESH_TYPE &m, MIDPOINT mid, typename MESH_TYPE::ScalarType thr = 0, bool RefineSelected = false)
        {
            EdgeLen <MESH_TYPE, typename MESH_TYPE::ScalarType> ep(thr);
            return RefineE(m, mid, ep, RefineSelected);
        }

    } // namespace tri
} // namespace tessellation
