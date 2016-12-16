#include "Quadtree.h"

#include <math.h>

void initQuadtree(Quadtree *qt, int height, int nbParticlesMin, int nbParticlesMax, 
				  double mMin, double mMax, double xMin, double xMax, double yMin, double yMax)
{
	qt->height = height;

	qt->nbCells = powl(4, height-1);
	qt->cells = (Cell *) malloc(qt->nbCells * sizeof(Particle));
	
	long dim = powl(2, height-1);
	double dX = (xMax - xMin) / dim;
	double dY = (yMax - yMin) / dim;

	double xMinCell = 0, yMinCell = 0;
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			initCell(qt->cells + (j * dim + i), nbParticlesMin / qt->nbCells, nbParticlesMax / qt->nbCells,
					 mMin, mMax, xMinCell, xMinCell+dX, yMinCell, yMinCell+dY);
			yMinCell += dY;
		}
		xMinCell += dX;
	}

	qt->centers = (Particle *) malloc(((powl(4, height-1) - 1) / 3) * sizeof(Particle));

}