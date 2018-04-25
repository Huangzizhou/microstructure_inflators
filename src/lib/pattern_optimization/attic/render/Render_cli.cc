// Use Mesa includes--should work on both on linux and Macs (with mesa installed // through macports)
#include <GL/osmesa.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>

#include <string> // sketchy png++ needs this... (but so do we)
#include <png++/png.hpp>
#include <stdexcept>

#include <boost/program_options.hpp>

#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldParser.hh>
#include <MeshFEM/EdgeFields.hh>
#include <LinearFEM.hh>

#include <CSGFEM/colors.hh>
#include <CSGFEM/draw.hh>

namespace po = boost::program_options;
using namespace std;

void usage(int exitVal, const po::options_description &visible_opts)
{
    cout << "Usage: Render_cli out.png mesh [edge_fields] [options]" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("outPng",      po::value<string>(), "output .png file")
        ("mesh",        po::value<string>(), "mesh to visualize")
        ("edgeFields",  po::value<string>(), "edge fields visualize")
        ;

    po::positional_options_description p;
    p.add("outPng",     1);
    p.add("mesh",       1);
    p.add("edgeFields", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help", "Produce this help message")
        ("edgeVectorField,v", po::value<string>(),        "Edge vector field")
        ("edgeScalarField,s", po::value<string>(),        "Edge scalar field")
        ("width,w",  po::value<int>()->default_value(1024), "output image width")
        ("height,h", po::value<int>()->default_value(768),  "output image height")
        ("viewBox", po::value<string>(),                  "view bounding box (Imagemagick format: WxH+X+Y)")
        ("key",                                           "draw colormap key")
        ("antialias,A",                                   "draw colormap key")
        ("wireframeWidth,W", po::value<double>()->default_value(0.0), "line width for wireframe")
        ("boundaryWidth,b",  po::value<double>()->default_value(1.0),  "line width for shaded boundary")
        ("arrowLineWidth,a", po::value<double>()->default_value(1.0),  "line width for vector field arrows")
        ("fieldScale,S",     po::value<double>(), "scale factor for scalar/vector fields (instead of normalizing)")
        ;

    po::options_description cli_opts;
    cli_opts.add(visible_opts).add(hidden_opts);
    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
    }

    if (vm.count("help"))
        usage(0, visible_opts);

    if ((vm.count("outPng") == 0) || (vm.count("mesh") == 0)) {
        cout << "Error: must specify output png and mesh files" << endl;
        usage(1, visible_opts);
    }

    return vm;
}

template<class _Mesh>
void drawMesh(const _Mesh &mesh) {
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < mesh.numElements(); ++i) {
        auto e = mesh.element(i);
        glVertex2f(e.vertex(0)->p[0], e.vertex(0)->p[1]);
        glVertex2f(e.vertex(1)->p[0], e.vertex(1)->p[1]);
        glVertex2f(e.vertex(2)->p[0], e.vertex(2)->p[1]);
    }
    glEnd();
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
    po::variables_map args = parseCmdLine(argc, argv);

    int width = args["width"].as<int>(), height = args["height"].as<int>();

    vector<MeshIO::IOVertex>  inVertices;
    vector<MeshIO::IOElement> inElements;
    string meshPath = args["mesh"].as<string>();
    auto type = load(meshPath, inVertices, inElements, MeshIO::FMT_GUESS,
                     MeshIO::MESH_GUESS);
    if (type != MeshIO::MESH_TRI)
        throw std::runtime_error("Only triangle meshes are supported");

    LinearFEM2D::Mesh<> mesh(inElements, inVertices);

    MSHFieldParser<2> *mshFields = NULL;
    if (MeshIO::guessFormat(meshPath) == MeshIO::FMT_MSH)
        mshFields = new MSHFieldParser<2>(meshPath);

    size_t numBE = mesh.numBoundaryElements();
    EdgeFields *edgeFields = NULL;
    vector<size_t> beToEdgeFieldIndex;
    if (args.count("edgeFields")) {
        edgeFields = new EdgeFields(args["edgeFields"].as<string>());
        auto incompatible = runtime_error("edgeFields doesn't match mesh");
        if (edgeFields->numEdges() != numBE) throw incompatible;
        for (size_t i = 0; i < mesh.numBoundaryElements(); ++i) {
            auto be = mesh.boundaryElement(i);
            beToEdgeFieldIndex.push_back(
                    edgeFields->edgeIndex(be.vertex(0).volumeVertex().index(),
                        be.vertex(1).volumeVertex().index()));
        }
    }

    OSMesaContext ctx = OSMesaCreateContextExt(OSMESA_RGBA, 16, 0, 0, NULL);
    unsigned char *osmesaBuffer = new unsigned char[4 * width * height];
    OSMesaMakeCurrent(ctx, osmesaBuffer, GL_UNSIGNED_BYTE, width, height);

    ////////////////////////////////////////////////////////////////////////////
    // Fit mesh in view (or take custom view from commandline)
    ////////////////////////////////////////////////////////////////////////////
    BBox<Vector2D> meshBox, frame;
    if (args.count("viewBox")) {
        throw std::runtime_error("Unimplemented");
    }
    else {
        meshBox = mesh.boundingBox();
    }

    double meshAspectRatio = meshBox.dimensions()[0] / meshBox.dimensions()[1];
    double viewAspectRatio = ((double) width) / height;
    frame = meshBox;
    // Expand the frame so that it has the same aspect ratio as the view.
    if (meshAspectRatio < viewAspectRatio)
         frame.expand(Vector2D(viewAspectRatio / meshAspectRatio - 1, 0));
    else frame.expand(Vector2D(0, meshAspectRatio / viewAspectRatio - 1));
    frame.expand(Vector2D(0.01, 0.01));

    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(frame.minCorner[0], frame.maxCorner[0],
            frame.minCorner[1], frame.maxCorner[1], -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    if (args.count("antialias")) {
        // Use antialiasing
        glEnable(GL_BLEND);
        // glBlendEquation(GL_FUNC_ADD);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_LINE_SMOOTH);
        glEnable(GL_POINT_SMOOTH);
    }

    // Draw mesh
    glColor3f(0.90, 0.90, 0.90);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    drawMesh(mesh);

    // Draw wireframe, if requested
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    if (args["wireframeWidth"].as<double>() > 0) {
        glLineWidth(args["wireframeWidth"].as<double>());
        glColor3f(0, 0, 0);
        drawMesh(mesh);
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
    // Draw boundary edges, potentially shaded with a scalar field.
    ScalarField<Real> edgeSField;
    ColorMap<RGBColorf, Real> colorMap(COLORMAP_JET);
    if (args.count("edgeScalarField")) {
        if (!edgeFields) throw std::runtime_error("No edge fields");
        edgeSField = edgeFields->field(args["edgeScalarField"].as<string>());
        colorMap.setRange(-std::abs(edgeSField.maxMag()),
                std::abs(edgeSField.maxMag()));
        glLineWidth(args["boundaryWidth"].as<double>());
    }
    glBegin(GL_LINES);
    for (size_t bei = 0; bei < numBE; ++bei) {
        if (edgeSField.domainSize() == numBE)
            glColor4fv(colorMap(edgeSField[beToEdgeFieldIndex[bei]]));
        else
            glColor3f(0.0, 0.0, 0.0);
        auto p0 = mesh.boundaryElement(bei).vertex(0).volumeVertex()->p;
        auto p1 = mesh.boundaryElement(bei).vertex(1).volumeVertex()->p;
        glVertex2f(p0[0], p0[1]);
        glVertex2f(p1[0], p1[1]);
    }
    glEnd();

    // Draw boundary vector field
    if (args.count("edgeVectorField")) {
        if (!edgeFields) throw std::runtime_error("No edge fields");
        VectorField<Real, 2> edgeVField = edgeFields->field(args["edgeVectorField"].as<string>());
        glLineWidth(args["arrowLineWidth"].as<double>());

        if (args.count("fieldScale") == 0) {
            // Normalize vector field for display
            Real maxNorm = edgeVField.maxMag();
            edgeVField *= 0.05 * meshBox.dimensions()[0] / maxNorm;
        }
        else {
            edgeVField *= args["fieldScale"].as<double>();
        }
        
        glColor3f(0.0, 0.0, 0.0);
        for (size_t bei = 0; bei < numBE; ++bei) {
            auto p0 = mesh.boundaryElement(bei).vertex(0).volumeVertex()->p;
            auto p1 = mesh.boundaryElement(bei).vertex(1).volumeVertex()->p;
            Vector2D midpoint = 0.5 * (p0 + p1);
            drawArrow2D(midpoint[0], midpoint[1],
                        edgeVField(bei)[0], edgeVField(bei)[1]);
        }
    }

    glFinish();

    png::image<png::rgba_pixel> pngWriter(width, height);
    for (size_t x = 0; x < width; ++x) {
        for (size_t y = 0; y < height; ++y) {
            pngWriter.set_pixel(x, y,
                    png::rgba_pixel(osmesaBuffer[4 * (width * (height - y - 1) + x) + 0],
                                    osmesaBuffer[4 * (width * (height - y - 1) + x) + 1],
                                    osmesaBuffer[4 * (width * (height - y - 1) + x) + 2],
                                    osmesaBuffer[4 * (width * (height - y - 1) + x) + 3]));
        }
    }

    pngWriter.write(args["outPng"].as<string>());

    OSMesaDestroyContext(ctx);
    delete osmesaBuffer;

    return 0;
}
