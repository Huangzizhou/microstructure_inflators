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

#include "glwidget.h"
#include "edge_mesh_type.h"
#include "tri_mesh_type.h"
#include <wrap/qt/trackball.h>
#include <wrap/gl/trimesh.h>

TwBar *bar,*sweepbar;
char * filename;          // filename of the mesh to load
vcg::Trackball track;     // the active manipulator

vcg::GlTrimesh<MyTriMesh> GLMesh;

MyTriMesh GuidanceM;

CMesh WireM;
CMesh SupportM;

std::string OutputMesh;

std::string WirePath="./flappy/bird.wire";//"./data/pattern1151_3x3x3.wire";
std::string DomainPath="./flappy/bird_guide.obj";//"./data/pattern1151_3x3x3_guide_mesh.obj";

float dir[3] = {-1, 1, -1};


void MoveToGround()
{
    vcg::tri::UpdateBounding<CMesh>::Box(WireM);
    for (size_t i=0;i<WireM.vert.size();i++)
        WireM.vert[i].P()-=WireM.bbox.min;

    for (size_t i=0;i<GuidanceM.vert.size();i++)
        GuidanceM.vert[i].P()-=WireM.bbox.min;

    vcg::tri::UpdateBounding<CMesh>::Box(WireM);
    vcg::tri::UpdateBounding<MyTriMesh>::Box(GuidanceM);
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

    vcg::tri::UpdateBounding<CMesh>::Box(WireM);
    vcg::tri::UpdateNormal<MyTriMesh>::PerFaceNormalized(GuidanceM);

    MoveToGround();
    WireM.printable=WireM.PrintabilityTest();
}

void LoadData()
{
    WireM.Load(WirePath);
    printf("Loaded Wire mesh with %d Vertices and %d Edges",WireM.vn,WireM.en);
    fflush(stdout);

    GuidanceM.Load(DomainPath);
    printf("Loaded Wire mesh with %d Vertices and %d Edges",WireM.vn,WireM.en);
    fflush(stdout);

    MoveToGround();

}

void TW_CALL SavePrintabilityData(void *)
{
    vcg::tri::io::ExporterOBJ<CMesh>::Save(WireM,"./out/wire_printable.wire",0);
    vcg::tri::io::ExporterOBJ<CMesh>::Save(SupportM,"./out/support.wire",0);
    vcg::tri::io::ExporterOBJ<MyTriMesh>::Save(GuidanceM,"./out/GuidanceM_transformed.obj",0);

}

void TW_CALL Rotate(void *)
{
    RotatePattern();
}

void TW_CALL InitPrintability(void *)
{
    vcg::tri::UpdateFlags<CMesh>::VertexClearS(WireM);
    WireM.SelectBoundaryVertices<MyTriMesh>(GuidanceM);
    WireM.printable=WireM.PrintabilityTest();
}

void TW_CALL AddPrintingSupport(void *)
{
    WireM.CreatePrintingSupport(SupportM,GuidanceM);
}

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{

    if (WireM.vn==0)
        LoadData();

    GLMesh.m=&GuidanceM;

    filename=0;
    //hasToPick=false;
    setWindowTitle(tr("Printability Optimizer GL"));
    bar = TwNewBar("Tools");
    TwDefine("Tools size='500 600' ");
    TwDefine("Tools position='40 40' ");
    TwAddButton(bar,"Rotate ",Rotate,0," label='Rotate'");

     // yellow
    // ...
    TwAddVarRW(bar, "RotDir", TW_TYPE_DIR3F, &dir, "");

    TwAddButton(bar,"Init Printability",InitPrintability,0," label='Init Printability'");
    TwAddButton(bar,"Add Support",AddPrintingSupport,0," label='Add Printing Support'");
    TwAddButton(bar,"Save",SavePrintabilityData,0,	" label='Save Data'");
}

void GLWidget::initializeGL ()
{
    glewInit();
    glClearColor(0, 0, 0, 0);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
}

void GLWidget::resizeGL (int w, int h)
{
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    TwWindowSize(w, h);
    initializeGL();
}

void GLWidget::paintGL ()
{
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, GLWidget::width()/(float)GLWidget::height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3.5f,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 1;
    track.GetView();
    glPushMatrix();
    track.Apply();
    glPushMatrix();

    vcg::glScale(3.0/WireM.bbox.Diag());
    glTranslate(-WireM.bbox.Center());
    WireM.GlDraw(0.3);
    SupportM.GlDraw(0.1,false);
    vcg::glColor(vcg::Color4b(0,0,0,255));
    GLMesh.Draw<vcg::GLW::DMWire,vcg::GLW::CMNone,vcg::GLW::TMNone>();

    glPopMatrix();
    //track.DrawPostApply();
    glPopMatrix();

    TwDraw();
}

void GLWidget::keyReleaseEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt) track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
    updateGL ();
}


void GLWidget::keyPressEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control) track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));

    TwKeyPressQt(e);
    updateGL ();
}

void GLWidget::mousePressEvent (QMouseEvent * e)
{
    if(!TwMousePressQt(this,e))
    {
        e->accept ();
        setFocus ();
        track.MouseDown(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG (e->button (), e->modifiers ()));
    }
    updateGL ();
}

void GLWidget::mouseMoveEvent (QMouseEvent * e)
{
    if (e->buttons ()) {
        track.MouseMove(QT2VCG_X(this, e), QT2VCG_Y(this, e));
        updateGL ();
    }
    TwMouseMotion(QTLogicalToDevice(this, e->x()), QTLogicalToDevice(this, e->y()));
}

void GLWidget::mouseReleaseEvent (QMouseEvent * e)
{
    track.MouseUp(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG(e->button (), e->modifiers ()));
    TwMouseReleaseQt(this,e);
    updateGL ();
}

void GLWidget::wheelEvent (QWheelEvent * e)
{
    const int WHEEL_STEP = 120;
    track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
    updateGL ();
}
