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

#include "edge_mesh_type.h"
#include <wrap/io_trimesh/export.h>
#include "tri_mesh_type.h"

std::string WirePath;
std::string GuidancePath;
std::string InflatedPath;
float dir[3];
float offset;
std::string OutputWire;
std::string OutputGuidance;
std::string OutputSupport;
std::string OutputInflated;

MyTriMesh GuidanceM;
MyTriMesh InflatedM;

CMesh WireM;
CMesh SupportM;

void InitPrintability()
{
    vcg::tri::UpdateFlags<CMesh>::VertexClearS(WireM);
    WireM.SelectBoundaryVertices<MyTriMesh>(GuidanceM);
    WireM.printable=WireM.PrintabilityTest();
}

void AddPrintingSupport()
{
    WireM.CreatePrintingSupport(SupportM,GuidanceM);
}

void MoveToGround(float bottom=0)
{
    vcg::tri::UpdateBounding<CMesh>::Box(WireM);
    for (size_t i=0;i<WireM.vert.size();i++)
        WireM.vert[i].P()-=WireM.bbox.min;

    for (size_t i=0;i<GuidanceM.vert.size();i++)
        GuidanceM.vert[i].P()-=WireM.bbox.min;

    for (size_t i=0;i<InflatedM.vert.size();i++)
        InflatedM.vert[i].P()-=WireM.bbox.min;

    for (size_t i=0;i<WireM.vert.size();i++)
        WireM.vert[i].P().Y()+=bottom;

    for (size_t i=0;i<GuidanceM.vert.size();i++)
        GuidanceM.vert[i].P().Y()+=bottom;

    for (size_t i=0;i<InflatedM.vert.size();i++)
        InflatedM.vert[i].P().Y()+=bottom;

    vcg::tri::UpdateBounding<CMesh>::Box(WireM);
    vcg::tri::UpdateBounding<MyTriMesh>::Box(GuidanceM);
    vcg::tri::UpdateBounding<MyTriMesh>::Box(InflatedM);
}

void RotatePattern()
{
    CoordType YDir(dir[0],dir[1],dir[2]);
    YDir.Normalize();
    vcg::Matrix33<ScalarType> Rot=vcg::RotationMatrix(CoordType(0,1,0),YDir);

    for (size_t i=0;i<GuidanceM.vert.size();i++)
        GuidanceM.vert[i].P()=Rot*GuidanceM.vert[i].P();

    for (size_t i=0;i<WireM.vert.size();i++)
        WireM.vert[i].P()=Rot*WireM.vert[i].P();

    for (size_t i=0;i<InflatedM.vert.size();i++)
        InflatedM.vert[i].P()=Rot*InflatedM.vert[i].P();

    vcg::tri::UpdateBounding<CMesh>::Box(WireM);
    vcg::tri::UpdateNormal<MyTriMesh>::PerFaceNormalized(GuidanceM);
    vcg::tri::UpdateNormal<MyTriMesh>::PerFaceNormalized(InflatedM);

    MoveToGround();
    //WireM.printable=WireM.PrintabilityTest();
}

int main(int argc, char *argv[])
{

    if (argc>1)
    {
        assert(argc==12);
        WirePath=std::string(argv[1]);
        GuidancePath=std::string(argv[2]);
        InflatedPath=std::string(argv[3]);

        dir[0]=atof(argv[4]);
        dir[1]=atof(argv[5]);
        dir[2]=atof(argv[6]);
        offset=atof(argv[7]);

        OutputWire=std::string(argv[8]);
        OutputSupport=std::string(argv[9]);
        OutputGuidance=std::string(argv[10]);
        OutputInflated=std::string(argv[11]);

        printf("Offesetting %5.5f \n",offset);
        fflush(stdout);

        WireM.Load(WirePath);
        printf("Loading %s as Input Wire Mesh\n",WirePath.c_str());

        GuidanceM.Load(GuidancePath);
        printf("Loading %s as Input Guidance Mesh\n",GuidancePath.c_str());
        fflush(stdout);

        InflatedM.Load(InflatedPath);
        printf("Loading %s as Inflated Mesh\n",InflatedPath.c_str());
        fflush(stdout);

        printf("Rotating along direction %5.5f %5.5f %5.5f\n",dir[0],dir[1],dir[2]);
        fflush(stdout);

        MoveToGround();

        RotatePattern();

        MoveToGround(offset);

        InitPrintability();
        if (!WireM.printable)
            printf("\nWARNING THE MESH IS NOT PRINTABLE\n");
        else
            printf("\nTHE MESH IS PRINTABLE\n");
        fflush(stdout);

        AddPrintingSupport();

        //MoveToGround();

        vcg::tri::io::ExporterOBJ<CMesh>::Save(WireM,OutputWire.c_str(),0);
        printf("Saving Output Wire Mesh in %s\n",OutputWire.c_str());

        vcg::tri::io::ExporterOBJ<CMesh>::Save(SupportM,OutputSupport.c_str(),0);
        printf("Saving Output Support Mesh in %s\n",OutputSupport.c_str());

        vcg::tri::io::ExporterOBJ<MyTriMesh>::Save(GuidanceM,OutputGuidance.c_str(),0);
        printf("Saving Output Guidance Mesh in %s\n",OutputGuidance.c_str());

        vcg::tri::io::ExporterOBJ<MyTriMesh>::Save(InflatedM,OutputInflated.c_str(),0);
        printf("Saving Output Inflated Mesh in %s\n",OutputInflated.c_str());

        fflush(stdout);
    }

}
