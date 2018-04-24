#include "WireInflator2D.h"
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <MeshFEM/MeshIO.hh>

#include <string>
#include <sstream>
#include <vector>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// James messed with the parameter range (WireMesh2D.h:setup),
// making the radius upper bound way too high.
// Manually choose a more reasonable one here.
double minRadius = 0.01, maxRadius = 0.1;
std::vector<std::pair<double, double>> range;


CellParameters rand_param(const WireInflator2D &wi)
{
	CellParameters p(wi.numberOfParameters());
	for (int i=0; i<int(wi.numberOfParameters()); ++i)
		p.parameter(i) = range[i].first + (double(rand()) / RAND_MAX) * (range[i].second - range[i].first);
	return p;
}

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cerr << "Usage:" << std::endl;
		std::cerr << argv[0] << " <wire_file> [<quads_file>]" << std::endl;
		return -1;
	}
	static const std::string wireMeshPath(argv[1]);

    auto wi = WireInflator2D::construct(wireMeshPath);
	WireInflator2D::OutMeshType mesh;

    CellParameters p_params = wi->createParameters();
    auto ops = wi->getParameterOperations();
    for (size_t i = 0; i < p_params.numberOfParameters(); ++i) {
        if (ops[i].type == ParameterOperation::Radius)
            range.push_back({minRadius, maxRadius});
        else range.push_back(wi->getParameterRange(i));
    }
    
    std::cout << p_params.numberOfParameters() << " parameters" << std::endl;

	// single cell generations
    if (0)
	{
		TessellationParameters t_params;

		for (size_t i=0; i<p_params.numberOfParameters(); ++i)
		{
			p_params.parameter(int(i)) = (range[i].first + range[i].second) / 2.0;
		}
		assert(wi->parametersValid(p_params));

		wi->generatePattern(p_params, t_params, mesh);

        std::vector<MeshIO::IOVertex> outVertices;
        std::vector<MeshIO::IOElement> outElements;
        for (const auto &p : mesh.nodes)
            outVertices.push_back(MeshIO::IOVertex(p[0], p[1], 0));
        for (const auto &e : mesh.elements)
            outElements.push_back(MeshIO::IOElement(e[0], e[1], e[2]));
        save("single.msh", outVertices, outElements);
	}

	// tiled grid generation
    if (0)
	{
		TessellationParameters t_params;

		const int w=10, h=10;
		Array2D<CellParameters *> grid(w,h);
		srand(time(NULL));
		for (int x=0; x<w; ++x)
			for (int y=0; y<h; ++y)
			{
				grid(x,y) = ((double(rand()) / RAND_MAX) < 0.2) ? NULL : new CellParameters();
				if (grid(x,y) != NULL)
				{
					*grid(x,y) = rand_param(*wi);
				}
			}

		wi->generateTiledPattern(grid, t_params, mesh);

		for (int x=0; x<10; ++x)
			for (int y=0; y<10; ++y)
				delete grid(x,y);

        std::vector<MeshIO::IOVertex> outVertices;
        std::vector<MeshIO::IOElement> outElements;
        for (const auto &p : mesh.nodes)
            outVertices.push_back(MeshIO::IOVertex(p[0], p[1], 0));
        for (const auto &e : mesh.elements)
            outElements.push_back(MeshIO::IOElement(e[0], e[1], e[2]));
        save("inflated.msh", outVertices, outElements);
	}

	if (argc <= 2)
		return 0;

	// quads generation
    // usage: pattern.obj quad_mesh.obj patterns.txt
    // patterns.txt is a subset of lines from the pattern lookup table
    // (with as many lines as quads in quad_mesh.obj).
    assert(argc == 4);
    // parse patterns.txt
    ifstream patternsFile(argv[3]);
    if (!patternsFile.is_open())
        throw std::runtime_error("Couldn't open " + string(argv[3]));
    std::string line;

    std::vector<CellParameters> quadParams;
    size_t nParams = p_params.numberOfParameters();
    while (std::getline(patternsFile, line)) {
        auto tokens = split(line, '\t');
        assert(tokens.size() == 1 + 3 + 1 + nParams);
        quadParams.emplace_back(wi->numberOfParameters());
        auto &p = quadParams.back();
        for (size_t i = 0; i < nParams; ++i)
            p.parameter(i) = std::stod(tokens[5 + i]);
    }

    // WARNING: CLIPPER MUST BE COMPILED WITH use_lines DEFINED FOR THIS!
    TessellationParameters t_params;
    static const bool averageThicknessOnBoundary = true;

    wi->generateQuadsPattern(std::string(argv[2]), quadParams, t_params,
                             mesh, averageThicknessOnBoundary);

    std::vector<MeshIO::IOVertex> outVertices;
    std::vector<MeshIO::IOElement> outElements;
    for (const auto &p : mesh.nodes)
        outVertices.push_back(MeshIO::IOVertex(p[0], p[1], 0));
    for (const auto &e : mesh.elements)
        outElements.push_back(MeshIO::IOElement(e[0], e[1], e[2]));
    save("tiled.msh", outVertices, outElements);

	return 0;
}
