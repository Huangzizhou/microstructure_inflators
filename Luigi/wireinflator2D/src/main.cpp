#include "WireInflator2D.h"
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <iostream>

#if 0
CellParameters rand_param(const PatternGen & pg)
{
	CellParameters p(pg.numberOfParameters());
	for (int i=0; i<int(pg.numberOfParameters()); ++i)
	{
		std::pair<double, double> range = pg.getParameterRange(i);
		p.parameter(i) = range.first + (double(rand()) / RAND_MAX) * (range.second - range.first);
	}
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

	// single cell generations
	{
		TessellationParameters t_params;
		CellParameters         p_params = wi->createParameters();
		std::cout << p_params.numberOfParameters() << " parameters" << std::endl;

		for (size_t i=0; i<p_params.numberOfParameters(); ++i)
		{
			std::pair<double,double> range = wi->patternGenerator().getParameterRange(int(i));
			p_params.parameter(int(i)) = (range.first + range.second) / 2.0;
		}
		assert(wi->patternGenerator().parametersValid(p_params));

		wi->generatePattern(p_params, t_params, mesh);
	}

	// tiled grid generation
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
					*grid(x,y) = rand_param(wi->patternGenerator());
				}
			}

		wi->generateTiledPattern(grid, t_params, mesh);

		for (int x=0; x<10; ++x)
			for (int y=0; y<10; ++y)
				delete grid(x,y);
	}

	if (argc <= 2)
		return 0;

	// quads generation
	static const std::string quadMeshFile(argv[2]);
	{
		TessellationParameters t_params;
		static const bool averageThicknessOnBoundary = true;

		srand(time(NULL));
		std::vector<CellParameters> quadParams;
		for (size_t i=0; i<15; ++i)
		{
			quadParams.push_back(rand_param(wi->patternGenerator()));
		}

		wi->generateQuadsPattern(quadMeshFile, quadParams, t_params, mesh, averageThicknessOnBoundary);
	}
	return 0;
}
#else

int main() {
    std::cout << "here" << std::endl;
    return 0;
}
#endif
