#ifndef QUADTREE_H
#define QUADTREE_H

#include "Particle.h"
#include "Cell.h"

struct Quadtree
{
	int height;
	int nbCells; // 4^(height-1)
	Cell *cells; // outer vertices of the tree, particle cells 
	int nbCenters; // (4^(height-1)-1)/3
	Particle *centers; // inner vertices of the tree, centers of mass
}
typedef Quadtree;

void initQuadtree(Quadtree *qt, int height, int nbParticlesMin, int nbParticlesMax, 
				  double mMin, double mMax, double xMin, double xMax, double yMin, double yMax);

#endif 