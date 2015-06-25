////////////////////////////////////////////////////////////////////////////////
// main.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      GLUT/AntTweakBar visualization of the parameterization
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//
//  Created:  12/04/2010 05:12:38
//  Revision History:
//      12/04/2010  Julian Panetta    Initial Revision
////////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <vector>
#include <string>
#include <limits>
#include <iostream>
#include <fstream>
#include <memory>
#include <AntTweakBar.h>

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#include <MeshIO.hh>
#include "Inflator.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Globals
////////////////////////////////////////////////////////////////////////////////
shared_ptr<ConstrainedInflator<2>> inflator;
vector<double> params;
float zoom;

////////////////////////////////////////////////////////////////////////////////
// IO Callback Routines
////////////////////////////////////////////////////////////////////////////////
void KeyboardFunc(unsigned char c, int x, int y) {
    TwEventKeyboardGLUT(c, x, y);   
    glutPostRedisplay();
}

void SpecialKeyboardFunc(int k, int x, int y) {
    TwEventSpecialGLUT(k, x, y);
    glutPostRedisplay();
}

void PassiveMotionFunc(int x, int y) {
    TwEventMouseMotionGLUT(x, y);
    glutPostRedisplay();
}

void MotionFunc(int x, int y) {
    TwEventMouseMotionGLUT(x, y);
    glutPostRedisplay();
}
    
void MouseFunc(int button, int state, int x, int y) {
    TwEventMouseButtonGLUT(button, state, x, y);
    glutPostRedisplay();
}

void inflate() {
    try {
        inflator->inflate(params);
    }
    catch (...) {
        inflator->clear();
        cerr << "Inflator fail." << endl;
    }
    glutPostRedisplay();
}

void TW_CALL setMaxVolCB(const void *value, void *info) {
    inflator->setMaxElementVolume(*(double *)value);
    inflate();
}

void TW_CALL getMaxVolCB(void *value, void *info) {
    *(double *)value = inflator->getMaxElementVolume();
}

void TW_CALL setParamCB(const void *value, void *idx) {
    params.at((size_t) idx) = *((const double *) value);
    inflate();
}

void TW_CALL getParamCB(void *value, void *idx) {
    *((double *) value) = params.at((size_t) idx);
}

void TW_CALL inflateCallback(void *caller) {
    inflate();
    glutPostRedisplay();
}

void TW_CALL saveCallback(void *caller) {
    inflate();
    MeshIO::save("out.msh", inflator->vertices(), inflator->elements());
    cout << "save!" << endl;
}

////////////////////////////////////////////////////////////////////////////////
// View setup
////////////////////////////////////////////////////////////////////////////////
void Reshape(int width, int height) {
    // // Set OpenGL viewport and camera
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // gluPerspective(40, ((double) width) / height, 1, 10);
    // glMatrixMode(GL_MODELVIEW);
    // glLoadIdentity();
    // gluLookAt(0,0,5, 0,0,0, 0,1,0);
    // glTranslatef(0, 0, -1);

    // Send the new window size to AntTweakBar
    TwWindowSize(width, height);

    glutPostRedisplay();
}

////////////////////////////////////////////////////////////////////////////////
/*! Applies the view transformation to the current matrix
*///////////////////////////////////////////////////////////////////////////////
void applyViewTransforms() {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

////////////////////////////////////////////////////////////////////////////////
/*! Called by GLUT when redisplay needed
*///////////////////////////////////////////////////////////////////////////////
void Display()
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    int width  = viewport[2];
    int height = viewport[3];

    glShadeModel(GL_FLAT);
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);

    // Use antialiasing
    glEnable(GL_BLEND);
    // glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);

    // Apply rotation/zoom
    applyViewTransforms();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glTranslatef(-0.5, -0.5, 0);

    glColor3f(1, 1, 1);
    glBegin(GL_TRIANGLES);
    for (const auto &e : inflator->elements()) {
        for (size_t i = 0; i < e.size(); ++i) {
            auto &v = inflator->vertices().at(e[i]);
            glVertex3f(v[0], v[1], v[2]);
        }
    }

    glEnd();

    TwDraw();
    glutSwapBuffers();
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on sucess)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    // Initialize GLUT
    glutInit(&argc, argv);

    int width = 1024, height = 768;
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutCreateWindow("Inflator");

    // Set GLUT event callbacks
    // - Directly redirect GLUT mouse button events to AntTweakBar
    glutMouseFunc(MouseFunc);
    // - Directly redirect GLUT mouse motion events to AntTweakBar
    glutMotionFunc(MotionFunc);
    // - Directly redirect GLUT mouse "passive" motion events to AntTweakBar
    glutPassiveMotionFunc(PassiveMotionFunc);
    // - Directly redirect GLUT key events to AntTweakBar
    glutKeyboardFunc(KeyboardFunc);
    // - Directly redirect GLUT special key events to AntTweakBar
    glutSpecialFunc(SpecialKeyboardFunc);
    // - Send 'glutGetModifers' function pointer to AntTweakBar;
    //   required because the GLUT key event functions do not report key
    //   modifiers states.
    TwGLUTModifiersFunc(glutGetModifiers);

    // Initialize AntTweakBar
    if (!TwInit(TW_OPENGL, NULL)) {
        // A fatal error occured    
        fprintf(stderr, "AntTweakBar initialization failed: %s\n",
                TwGetLastError());
        return 1;
    }

    TwBar *bar = TwNewBar("Parameters");
    TwDefine(" Parameters color='128 128 128' alpha='192' ");

    // Set GLUT callbacks
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);


    if (argc != 2) {
        cerr << "Usage: inflate_viewer wire_path" << endl;
        exit(-1);
    }

    string wirePath(argv[1]);

    vector<string> constraints;
    inflator = make_shared<ConstrainedInflator<2>>(constraints, wirePath);
    params.resize(inflator->numParameters());

    // Parameter specification precedence:
    //  --parameters presides over .dof, which presides over .opt, which
    //  presides over defaults.
    Real defaultOffset = 0.0;
    Real defaultThickness = 0.1;
    for (size_t p = 0; p < params.size(); ++p) {
        const char *rangeString;
        if (inflator->parameterType(p) == ParameterType::Thickness) {
            params[p] = defaultThickness;
            rangeString = " min=0.01 max=5 step=0.01";
        }
        else {
            params[p] = defaultOffset;
            rangeString = " min=-0.5 max=0.5 step=0.01";
        }

        TwAddVarCB(bar, ("Param " + to_string(p)).c_str(), TW_TYPE_DOUBLE,
                setParamCB, getParamCB, (void *) p, rangeString);
    }

    TwAddVarCB(bar, "Max Element Vol", TW_TYPE_DOUBLE,
            setMaxVolCB, getMaxVolCB, NULL, " min=0.0001 max=0.01 step=0.0001");

    TwAddButton(bar, "Inflate", inflateCallback, NULL, "key=i");
    TwAddButton(bar, "Save", saveCallback, NULL, "key=s");

    // Call the GLUT main loop
    glutMainLoop();

    return 0;
}
