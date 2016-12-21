#include "Cell.h"
#include "Particle.h"
#include "Quadtree.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#define EPS 1e-15

// TESTING HELPER METHODS

// Compute the maximum and mean relative errors on forces applied on a group of cells.
// compared with the forces computed naively in a unique global cell
// relative erros are d(fx)/fx, and d(fy)/fy
void computeRelativeErrors(Cell *cMerged, int nbCells, Cell *cells, double *maxRelErr, double *meanRelErr)
{
	double error = 0;
	*maxRelErr = 0;
	*meanRelErr = 0;
	Particle *p1, *p2;
	Cell *c;

	int shift = 0;
	for (int cellNo = 0; cellNo < nbCells; cellNo++)
	{
		c = cells + cellNo;
		for (int partNo = 0; partNo < c->nbParticles; partNo++)
		{
			p1 = c->particles + partNo;
			p2 = cMerged->particles + (shift + partNo);
			error = fabs(p1->fx - p2->fx) / fabs(p2->fx);
			*maxRelErr = (error > *maxRelErr) ? error : *maxRelErr;
			*meanRelErr += error;
			error = fabs(p1->fy - p2->fy) / fabs(p2->fy);
			*maxRelErr = (error > *maxRelErr) ? error : *maxRelErr;
			*meanRelErr += error;
		}
		shift += c->nbParticles;
	}

	*meanRelErr /= (2*cMerged->nbParticles);
}

// Check that the cell c has the specified bounding box (with eps precision)
void checkCell(Cell *c, double xMin, double xMax, double yMin, double yMax, double eps)
{
	assert(fabs(c->xMin - xMin) < eps && 
		   fabs(c->xMax - xMax) < eps &&
		   fabs(c->yMin - yMin) < eps &&
		   fabs(c->yMax - yMax) < eps); 
}




// TESTS

void testPonP()
{
	printf("Testing inter-particle interactions computation: ");

	Particle p1, p2;
	initParticle(&p1, 2.0, 3.0, 4.0);
	initParticle(&p2, 2.5, 3.5, 4.5);

	PonP(&p1, &p2);
	double expected = -4.7164e-10;
	assert(p1.fx == 0 && p1.fx == 0 && fabs(p2.fx - expected) < EPS && fabs(p2.fy - expected) < EPS);
	
	printf("OK\n");
}

void testP2M()
{
	printf("Testing center of mass computation: ");
	Cell c;
	c.nbParticles = 2;
	c.particles = (Particle *) malloc(c.nbParticles * sizeof(Particle));
	initParticle(c.particles, 2.0, 3.0, 4.0);
	initParticle(c.particles + 1, 2.5, 3.5, 4.5);

	Particle m;
	P2M(&m, &c);
	assert(fabs(m.m - 4.5) < EPS && fabs(m.x - 14.75/4.5) < EPS && fabs(m.y - 19.25/4.5) < EPS);

	free(c.particles);

	printf("OK\n");
}


void testM2P()
{
	printf("Testing M2P: ");
	
	Cell c[2];
	Cell cMerged;
	Particle cm0, cm1;

	srand(42);
	initCell(c+0, 50, 100, 1e30, 1e32, 0, 1e17, 0, 1e17);
	initCell(c+1, 50, 100, 1e30, 1e32, 0e17, 1e17, 3e17, 4e17);
	mergeCell(&cMerged, 2, c);
	P2M(&cm0, c+0);
	P2M(&cm1, c+1);

	P2P_in(&cMerged);

	M2P(&cm0, c+1);
	M2P(&cm1, c+0);
	P2P_in(c+0);
	P2P_in(c+1);

	double maxRelativeError;
	double meanRelativeError;
	computeRelativeErrors(&cMerged, 2, c, &maxRelativeError, &meanRelativeError);

	freeCell(c+0);
	freeCell(c+1);
	freeCell(&cMerged);

	printf("max relative error = %e, average relative error = %e \n", maxRelativeError, meanRelativeError);
}


void testInitQuadtree()
{
	printf("Testing quadtree initialization (w/ Morton index): ");
	
	srand(42);
	Quadtree qt;
	initQuadtree(&qt, 3, 10, 20, 1.0, 2.0, 0.0, 1.0, 0.0, 1.0);

	checkCell(qt.cells + 0, 0.0, 0.25, 0.0, 0.25, EPS);
	checkCell(qt.cells + 1, 0.25, 0.5, 0.0, 0.25, EPS);
	checkCell(qt.cells + 2, 0.0, 0.25, 0.25, 0.5, EPS);
	checkCell(qt.cells + 3, 0.25, 0.5, 0.25, 0.5, EPS);
	checkCell(qt.cells + 4, 0.5, 0.75, 0.0, 0.25, EPS);
	checkCell(qt.cells + 5, 0.75, 1.0, 0.0, 0.25, EPS);
	checkCell(qt.cells + 6, 0.5, 0.75, 0.25, 0.5, EPS);
	checkCell(qt.cells + 7, 0.75, 1.0, 0.25, 0.5, EPS);
	checkCell(qt.cells + 8, 0.0, 0.25, 0.5, 0.75, EPS);
	checkCell(qt.cells + 9, 0.25, 0.5, 0.5, 0.75, EPS);
	checkCell(qt.cells + 10, 0.0, 0.25, 0.75, 1.0, EPS);
	checkCell(qt.cells + 11, 0.25, 0.5, 0.75, 1.0, EPS);
	checkCell(qt.cells + 12, 0.5, 0.75, 0.5, 0.75, EPS);
	checkCell(qt.cells + 13, 0.75, 1.0, 0.5, 0.75, EPS);
	checkCell(qt.cells + 14, 0.5, 0.75, 0.75, 1.0, EPS);
	checkCell(qt.cells + 15, 0.75, 1.0, 0.75, 1.0, EPS);

	freeQuadtree(&qt);

	printf("OK\n");
}




void testRegression1()
{
	srand(42);
	printf("Regression test 1, comparing results of P2P_in, P2P_ext on a 4x4 grid VS "
		   "P2P_in on one global cell: ");

	Cell cMerged;
	Quadtree qt;
	initQuadtree(&qt, 3, 10000, 20000, 1e30, 1e32, 0, 1e17, 0, 1e17);
	mergeCell(&cMerged, qt.nbCells, qt.cells);

	P2P_in(&cMerged);
	
	for (int i = 0; i < qt.nbCells; i++)
		for (int j = 0; j < qt.nbCells; j++)
		{
			if (i == j)
				P2P_in(qt.cells + i);
			else
				P2P_ext(qt.cells+j, qt.cells+i);
		}

	double maxRelativeError;
	double meanRelativeError;
	computeRelativeErrors(&cMerged, qt.nbCells, qt.cells, &maxRelativeError, &meanRelativeError);

	printf("max relative error = %e, average relative error = %e \n", maxRelativeError, meanRelativeError);

	freeCell(&cMerged);
	freeQuadtree(&qt);
}

void testRegression2()
{
	srand(42);
	printf("Regression test 2, comparing results of a Barnes-Hut like method on a 4x4 grid VS "
		   "P2P_in on one global cell: ");

	Cell cMerged;
	Quadtree qt;
	initQuadtree(&qt, 3, 10000, 20000, 1e30, 1e32, 0, 1e17, 0, 1e17);
	mergeCell(&cMerged, qt.nbCells, qt.cells);

	P2P_in(&cMerged);

	computeCMs(&qt);
	computeAllForces(&qt);	

	double maxRelativeError;
	double meanRelativeError;
	computeRelativeErrors(&cMerged, qt.nbCells, qt.cells, &maxRelativeError, &meanRelativeError);

	printf("max relative error = %e, average relative error = %e \n", maxRelativeError, meanRelativeError);

	freeCell(&cMerged);
	freeQuadtree(&qt);
}



int main(int argc, char const *argv[])
{
	testPonP();
	testP2M();
	testInitQuadtree();
	testM2P();
	testRegression1();
	testRegression2();

	return EXIT_SUCCESS;
}