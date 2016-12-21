#include "Morton.h"
#include "Quadtree.h"

#include <math.h>

#define DIST(x1, y1, x2, y2) sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))
#define FAR_FIELD_LIMIT 0.707

// Private auxiliary functions prototypes

// Returns the distance between the cell c and the center of mass cm
double distance(Particle *cm, Cell *c);

// Compute the centers of mass of the i-th vertex and its children
void computeCMrec(Quadtree *qt, int cmNo);

// Compute the gravitationnal force exerted on each particle of the
// i-th cell by the j-th center of mass
// l is width of the area approximated by the cm
void computeForcesRec(Quadtree *qt, int cellNo, int cmNo, double l);


// Public functions

// Initialise the cells of a quadtree of specified height, built on the 2D area [xMin, xMax]*[yMin, yMax]
// containing a total of particles picked randomly in [[nbPartMin, nbPartMax]], 
// with masses picked randomly in [[mMin, mMax]]   
void initQuadtree(Quadtree *qt, int height, int nbPartMin, int nbPartMax, 
				  double mMin, double mMax, double xMin, double xMax, double yMin, double yMax)
{
	qt->height = height;
	qt->xMin = xMin;
	qt->xMax = xMax;
	qt->yMin = yMin;
	qt->yMax = yMax;
	qt->nbCells = powl(4, height-1);
	qt->cells = (Cell *) malloc(qt->nbCells * sizeof(Cell));
	qt->nbCMs = (powl(4, height) - 1) / 3;
	qt->CMs = (Particle *) malloc(qt->nbCMs * sizeof(Particle));
	qt->firstOuterCM = qt->nbCMs - qt->nbCells; 
	
	long dim = powl(2, height-1);
	double dX = (xMax - xMin) / (double)dim;
	double dY = (yMax - yMin) / (double)dim;
	int nbPartPerCellMin = nbPartMin / qt->nbCells;
	int nbPartPerCellMax = nbPartMax / qt->nbCells;

	int cellNo;
	for (int x = 0; x < dim; x++)
		for (int y = 0; y < dim; y++)
		{
			cellNo = xy_to_morton(x, y);
			initCell(qt->cells + cellNo, nbPartPerCellMin, nbPartPerCellMax, mMin, mMax, 
					 x*dX, (x+1)*dX, y*dY, (y+1)*dY);
		}
}

// Release the ressources associated with the specified quadtree
void freeQuadtree(Quadtree *qt)
{
	for (int i = 0; i < qt->nbCells; i++)
		freeCell(qt->cells + i);
	free(qt->cells);
	free(qt->CMs);
}


// Compute all the centers of mass of the specified quadtree recursively, 
// starting by the lower level ones
void computeCMs(Quadtree *qt)
{
	computeCMrec(qt, 0);
}

// Compute the gravitational force exerted on each particle of the i-th cell
void computeForces(Quadtree *qt, int cellNo)
{
	computeForcesRec(qt, cellNo, 0, qt->xMax - qt->xMin);
}

// Compute the gravitationnal force exerted on each particule of the quadtree
void computeAllForces(Quadtree *qt)
{
	for (int i = 0; i < qt->nbCells; i++)
		computeForces(qt, i);
}


// Private auxiliary functions

double distance(Particle *cm, Cell *c)
{
	double d = 0;

	if (cm->x < c->xMin)
	{
		if (cm->y < c->yMin)
			d = DIST(c->xMin, c->yMin, cm->x, cm->y);
		else if (cm->y > c->yMax)
			d = DIST(c->xMin, c->yMax, cm->x, cm->y);
		else
			d = c->xMin - cm->x;
	}
	else if (cm->x > c->xMax)
	{
		if (cm->y < c->yMin)
			d = DIST(c->xMax, c->yMin, cm->x, cm->y);
		else if (cm->y > c->yMax)
			d = DIST(c->xMax, c->yMax, cm->x, cm->y);
		else
			d = cm->x - c->xMax;
	}
	else
	{
		if (cm->y < c->yMin)
			d = c->yMin - cm->y;
		else if (cm->y > c->yMax)
			d = cm->y - c->yMax;
	}

	return d;
}

void computeCMrec(Quadtree *qt, int cmNo)
{
	if (cmNo < qt->firstOuterCM)
	{
		computeCMrec(qt, 4*cmNo+1);
		computeCMrec(qt, 4*cmNo+2);
		computeCMrec(qt, 4*cmNo+3);
		computeCMrec(qt, 4*cmNo+4);
		M2M(qt->CMs + cmNo, qt->CMs + 4*cmNo+1);
	}
	else
	{
		P2M(qt->CMs + cmNo, qt->cells + (cmNo - qt->firstOuterCM));
	}
}

void computeForcesRec(Quadtree *qt, int cellNo, int cmNo, double l)
{
	// We consider that a center of mass is in the far field of a cell
	// if d/w > 0.707
	// where d is the distance between the cm and the cell
	// and w is the width of the area approximated by the cm
	double d = distance(qt->CMs + cmNo, qt->cells + cellNo); 

	if (d > l / FAR_FIELD_LIMIT)
	{
		// printf("M2P from cm %d (width = %e, distance %e) on cell %d\n", cmNo, l, d, cellNo);
		M2P(qt->CMs + cmNo, qt->cells + cellNo);
	}
	else if (cmNo < qt->firstOuterCM)
	{
		computeForcesRec(qt, cellNo, 4*cmNo+1, l/2);
		computeForcesRec(qt, cellNo, 4*cmNo+2, l/2);
		computeForcesRec(qt, cellNo, 4*cmNo+3, l/2);
		computeForcesRec(qt, cellNo, 4*cmNo+4, l/2);
	}
	else if (cmNo == qt->firstOuterCM + cellNo)
	{
		// printf("P2P_in on cell %d\n", cellNo);
		P2P_in(qt->cells + cellNo);
	}
	else
	{
		// printf("P2P_ext from cell %d on cell %d\n", cmNo - qt->firstOuterCM, cellNo);
		P2P_ext(qt->cells + (cmNo - qt->firstOuterCM), qt->cells + cellNo);
	}
}

