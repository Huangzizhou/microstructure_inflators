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

#ifndef TRI_MESH_TYPE_H
#define TRI_MESH_TYPE_H

#include <vcg/space/intersection3.h>

//the guidance mesh
class MyTriFace;
class MyTriVertex;
struct TriUsedTypes: public vcg::UsedTypes<vcg::Use<MyTriVertex>::AsVertexType,
        vcg::Use<MyTriFace>::AsFaceType>{};


class MyTriVertex:public vcg::Vertex<TriUsedTypes,
        vcg::vertex::Coord3f,
        vcg::vertex::Normal3f,
        vcg::vertex::Mark,
        vcg::vertex::BitFlags>
{};

class MyTriFace:public vcg::Face<TriUsedTypes,
        vcg::face::VertexRef,
        vcg::face::FFAdj,
        vcg::face::BitFlags,
        vcg::face::Normal3f>
{};


class MyTriMesh: public vcg::tri::TriMesh< std::vector<MyTriVertex>,
        std::vector<MyTriFace > >
{
public:
    void Load(std::string &path)
    {
        int mask;
        vcg::tri::io::ImporterOBJ<MyTriMesh>::LoadMask(path.c_str(),mask);
        vcg::tri::io::ImporterOBJ<MyTriMesh>::Open(*this,path.c_str(),mask);
        vcg::tri::UpdateBounding<MyTriMesh>::Box(*this);
        vcg::tri::UpdateNormal<MyTriMesh>::PerFaceNormalized(*this);
    }

    bool Intersect(CoordType pos,
                   CoordType direction,
                   int &faceNum)
    {
        for (size_t i=0;i<face.size();i++)
        {
            vcg::Ray3<ScalarType> R3(pos,direction);
            CoordType p[3];
            p[0]=face[i].P(0);
            p[1]=face[i].P(1);
            p[2]=face[i].P(2);
            ScalarType t,alpha,beta,gamma;
            if (vcg::IntersectionRayTriangle(R3,p[0],p[1],p[2],t,alpha,beta))
            {
                gamma=1-(alpha+beta);
                CoordType NormTris=face[i].N();
                if ((NormTris*direction)<(-0.1))
                {
                    faceNum=i;
                    return true;
                }
            }
        }
        return false;
    }

};



#endif
