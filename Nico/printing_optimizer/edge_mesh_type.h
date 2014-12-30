/****************************************************************************
**
** Copyright (C) 2011 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef EDGE_MESH_TYPE_H
#define EDGE_MESH_TYPE_H

#ifdef USE_OPENGL
#include <GL/glew.h>
#include <QGLWidget>
#include <QKeyEvent>
#endif

//vcg imports
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <wrap/io_trimesh/import.h>

//wrapper imports
#ifdef USE_OPENGL
#include <wrap/io_trimesh/import.h>
//#include <wrap/gl/space.h>
//#include <wrap/gl/math.h>
#include <wrap/gl/trimesh.h>
#include <wrap/gui/trackball.h>
#include <wrap/gl/addons.h>
#include <wrap/io_trimesh/export.h>
#include <QDir>
#endif

#include <algorithm>
#include <vcg/space/tetra3.h>
#include <vcg/space/distance3.h>
#include <vcg/space/line3.h>

typedef float ScalarType;
typedef vcg::Point3<ScalarType> CoordType;

using namespace vcg;
class CEdge;
class CVertex;


struct MyUsedTypes : public UsedTypes<	Use<CVertex>::AsVertexType,Use<CEdge>::AsEdgeType>{};

// compositing wanted proprieties
class CVertex : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f,
        vcg::vertex::Normal3f, vcg::vertex::BitFlags,
        vcg::vertex::VEAdj,vcg::vertex::Color4b>
{
public:
};

class CEdge   : public vcg::Edge<MyUsedTypes, edge::VertexRef,
                                 edge::BitFlags,edge::VEAdj>
{
public:
};

class CMesh   : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CEdge> >
{
public:
    bool printable;

#ifdef USE_OPENGL
    void GlDraw(ScalarType size=0.02,bool scaleEdge=true)
    {

        for (size_t i=0;i<vert.size();i++)
        {
            vcg::Color4b currCol;

            if (vert[i].IsS())
                currCol=vcg::Color4b(255,0,0,255);
            else
                currCol=vcg::Color4b(0,0,255,255);

            vcg::glColor(currCol);
            vcg::Add_Ons::glPoint<vcg::Add_Ons::DMSolid>(vert[i].P(),size,4,4);
        }

        if (printable)
            vcg::glColor(vcg::Color4b(0,255,0,255));
        else
            vcg::glColor(vcg::Color4b(255,255,0,255));

        for (size_t i=0;i<edge.size();i++)
        {
            CoordType p0=edge[i].V(0)->P();
            CoordType p1=edge[i].V(1)->P();
            CoordType bary=(p0+p1)/2.0;
            if (scaleEdge)
            {
                p0=p0*0.9+bary*0.1;
                p1=p1*0.9+bary*0.1;
            }
            //CoordType dir=p1-p0;
            vcg::Add_Ons::glCylinder<vcg::Add_Ons::DMSolid>(p0,p1,size,4);
        }
    }

#endif

    template <class TriMesh>
    void SelectBoundaryVertices(TriMesh &trimesh)
    {
        ScalarType EPS=0.00001;
        for (size_t i=0;i<vert.size();i++)
            for (size_t j=0;j<trimesh.face.size();j++)
            {
                vcg::Triangle3<ScalarType> T3;

                T3.P(0)=trimesh.face[j].cP(0);
                T3.P(1)=trimesh.face[j].cP(1);
                T3.P(2)=trimesh.face[j].cP(2);

                CoordType clos;
                ScalarType dist;
                vcg::TrianglePointDistance(T3,vert[i].P(),dist,clos);
                if ((dist<EPS)&&(trimesh.face[j].N().Y()<-EPS))vert[i].SetS();
            }
    }

    template <class TriMesh>
    void CreatePrintingSupport(CMesh &Support,
                               TriMesh &Guidance)
    {
        Support.Clear();
        for (size_t i=0;i<vert.size();i++)
        {
            if (!vert[i].IsS())continue;
            int faceNum;
            if (Guidance.Intersect( vert[i].P(),CoordType(0,-1,0),faceNum))continue;

            CoordType pos0=vert[i].P();
            CoordType pos1=vert[i].P();
            pos1.Y()=0;

            vcg::tri::Allocator<CMesh>::AddVertex(Support,pos0);
            vcg::tri::Allocator<CMesh>::AddVertex(Support,pos1);
            vcg::tri::Allocator<CMesh>::AddEdges(Support,1);
            int size=Support.vert.size();
            Support.edge.back().V(0)=&Support.vert[size-1];
            Support.edge.back().V(1)=&Support.vert[size-2];
        }
        Support.printable=true;
    }

    bool PrintabilityTest()
    {
        ScalarType EPS=0.0001;
        //vcg::tri::UpdateBounding<CMesh>::Box(*this);
        std::vector<int> status(vert.size(),0); //1 printable 0 unknown

        //calculate the star of vertices
        std::vector<std::vector<std::pair<int,bool> > > SupportN(vert.size(),std::vector<std::pair<int,bool> >());
        for (size_t i=0;i<edge.size();i++)
        {
            int N0=vcg::tri::Index(*this,edge[i].V(0));
            int N1=vcg::tri::Index(*this,edge[i].V(1));
            //check the verse of the edge
            CoordType P0=vert[N0].P();
            CoordType P1=vert[N1].P();
            if (fabs(P0.Y()-P1.Y())<EPS)
            {
                SupportN[N0].push_back(std::pair<int,bool>(N1,false));
                SupportN[N1].push_back(std::pair<int,bool>(N0,false));
            }
            else
            if (P0.Y()<P1.Y())
                SupportN[N1].push_back(std::pair<int,bool>(N0,true));
            else
                SupportN[N0].push_back(std::pair<int,bool>(N1,true));
        }

        //then set initialize bottom vertices
        for (size_t i=0;i<vert.size();i++)
        {
            if (vert[i].IsS())
            {
               status[i]=1;
               continue;
            }
            //there's no way this node to have support
            if ((SupportN[i].size())==0) return false;

            if ((SupportN[i].size())>0)
            {
                for (size_t j=0;j<SupportN[i].size();j++)
                {
                    if (SupportN[i][j].second)
                    {
                        status[i]=1;
                        break;
                    }
                }
            }
        }

        //bool IsPrint=true;
        std::vector<int> StackPrint;
        for (size_t i=0;i<status.size();i++)
        {
            if(status[i]==1)StackPrint.push_back(i);
        }
        if (StackPrint.size()==status.size()) return true;

//        printf("BFS check\n");
//        fflush(stdout);
        //BFS
        do
        {
         int currentN=StackPrint.back();
         assert(status[currentN]==1);
         StackPrint.pop_back();
         for (size_t i=0;i<SupportN[currentN].size();i++)
         {
             int HorizontalN=SupportN[currentN][i].first;

             if (status[HorizontalN]==1)continue;
             if (SupportN[currentN][i].second)continue;//only horizontal
             status[HorizontalN]=1;
             StackPrint.push_back(HorizontalN);
         }
        }while (!StackPrint.empty());

        for (size_t i=0;i<status.size();i++)
         if (status[i]==0)return false;

        return true;
    }

    void Load(std::string &path)
    {
        int mask;
        vcg::tri::io::ImporterOBJ<CMesh>::LoadMask(path.c_str(),mask);
        vcg::tri::io::ImporterOBJ<CMesh>::Open(*this,path.c_str(),mask);
        vcg::tri::UpdateBounding<CMesh>::Box(*this);
    }
};
#endif
