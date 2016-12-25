#ifndef QUADTREE_H
#define QUADTREE_H

#include "Multipole.h"
#include "Cell.h"

// Perfect complete 4-ary tree
// Below each leaf is stored a Cell, square subdivision of space containing particles.
// Each leaf is the center of mass of the particles contained in the corresponding Cell.
// Each inner vertex is the center of mass of its children vertices. 
struct Quadtree
{
	// properties of the quadtree
	int height;
	double xMin;
	double xMax;
	double yMin;
	double yMax;

	// there is one cell per outer vertex of the tree : 4^(height-1)
	int nbCells; 
	Cell *cells;
	
	// there is one center of mass per vertex of the tree : (4^(height-1)-1)/3 
	int nbMultipoles;
	Multipole *multipoles;

	// Index of the first outer vertex of the tree : nbCMS - nbCells
	int firstOuterCM;  
}
typedef Quadtree;

// Initialise the cells of a quadtree of specified height, built on the 2D area [xMin, xMax]*[yMin, yMax]
// containing a total of particles picked randomly in [[nbPartMin, nbPartMax]], 
// with masses picked randomly in [[mMin, mMax]]   
void initQuadtree(Quadtree *qt, int height, int nbPartMin, int nbPartMax, 
				  double mMin, double mMax, double xMin, double xMax, double yMin, double yMax);

// Release the ressources associated with the specified quadtree
void freeQuadtree(Quadtree *qt);

// Compute all the centers of mass of the specified quadtree recursively, 
// starting by the lower level ones
void computeMultipoles(Quadtree *qt);

// Compute the gravitationnal force exerted on each particule of the quadtree
void computeForces(Quadtree *qt);

#endif 