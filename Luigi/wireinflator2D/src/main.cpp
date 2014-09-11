#include "WireInflator2D.h"
#include <cassert>

typedef WireInflator2D::PatternGen PatternGen;

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
	(void)argc;
	(void)argv;

	static const std::string wireMeshPath = "../meshes/octa_cell.obj";

	WireInflator2D wi(wireMeshPath);

	// single cell generations
	{
		TessellationParameters t_params;
		CellParameters         p_params = wi.createParameters();

		for (size_t i=0; i<p_params.numberOfParameters(); ++i)
		{
			std::pair<double,double> range = wi.patternGenerator().getParameterRange(int(i));
			p_params.parameter(int(i)) = (range.first + range.second) / 2.0;
		}
		assert(wi.patternGenerator().parametersValid(p_params));

		WireInflator2D::OutMeshType mesh;
		wi.generatePattern(p_params, t_params, mesh);
	}

	// tiled grid generation
	{
		TessellationParameters t_params;

		const int w=10, h=10;
		Array2D<CellParameters *> grid(w,h);
		for (int x=0; x<w; ++x)
			for (int y=0; y<h; ++y)
			{
				grid(x,y) = ((double(rand()) / RAND_MAX) < 0.2) ? NULL : new CellParameters();
				if (grid(x,y) != NULL)
				{
					*grid(x,y) = rand_param(wi.patternGenerator());
				}
			}

		WireInflator2D::OutMeshType mesh;
		wi.generateTiledPattern(grid, t_params, mesh);

		for (int x=0; x<10; ++x)
			for (int y=0; y<10; ++y)
				delete grid(x,y);
	}
	return 0;
}
