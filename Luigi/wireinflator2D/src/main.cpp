#include "WireInflator2D.h"
#include <cassert>

WireInflator2D::PatternParameters rand_param(void)
{
	typedef WireInflator2D::PatternParameters PatternParameters;

	PatternParameters p;
	for (int i=0; i<int(p.numberOfParameters()); ++i)
	{
		std::pair<double, double> range = p.parameterRange(i);
		p.parameter(i) = range.first + (double(rand()) / RAND_MAX) * (range.second - range.first);
	}
	return p;
}

int main(int argc, char* argv[])
{
	(void)argc;
	(void)argv;

	// single cell generations
	{
		TessellationParameters            t_params;
		WireInflator2D::PatternParameters p_params;
		for (size_t i=0; i<p_params.numberOfParameters(); ++i)
		{
			std::pair<double,double> range = p_params.parameterRange(int(i));
			p_params.parameter(int(i)) = (range.first + range.second) / 2.0;
		}
		assert(p_params.isValid());
		WireInflator2D::OutMeshType mesh;
		WireInflator2D::generatePattern(p_params, t_params, mesh);
	}

	// tiled grid generation
	{
		TessellationParameters t_params;

		size_t w=10, h=10;
		Array2D<WireInflator2D::PatternParameters *> grid(w,h);
		for (int x=0; x<10; ++x)
			for (int y=0; y<10; ++y)
			{
				grid(x,y) = ((double(rand()) / RAND_MAX) < 0.2) ? NULL : new WireInflator2D::PatternParameters();
				if (grid(x,y) != NULL)
				{
					*grid(x,y) = rand_param();
				}
			}

		WireInflator2D::OutMeshType mesh;
		WireInflator2D::generateTiledPattern(grid, t_params, mesh);

		for (int x=0; x<10; ++x)
			for (int y=0; y<10; ++y)
				delete grid(x,y);
	}
	return 0;
}
