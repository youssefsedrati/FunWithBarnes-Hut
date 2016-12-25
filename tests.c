#include "Cell.h"
#include "Quadtree.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


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
	Cell *c;

	int shift = 0;
	for (int cellNo = 0; cellNo < nbCells; cellNo++)
	{
		c = cells + cellNo;
		for (int partNo = 0; partNo < c->nbParticles; partNo++)
		{
			error = fabs(c->fx[partNo] - cMerged->fx[shift + partNo]) / fabs(cMerged->fx[shift + partNo]);
			*maxRelErr = (error > *maxRelErr) ? error : *maxRelErr;
			*meanRelErr += error;
			error = fabs(c->fy[partNo] - cMerged->fy[shift + partNo]) / fabs(cMerged->fy[shift + partNo]);
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


// UNIT TESTS

void testP2M()
{
	printf("Testing center of mass computation: ");
	Cell c;
	initCell(&c, 2, 2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
	c.nbParticles = 2;
	c.m[0] = 2.0;
	c.m[1] = 2.5;
	c.x[0] = 3.0;
	c.x[1] = 3.5;
	c.y[0] = 4.0;
	c.y[1] = 4.5;
	
	Multipole m;
	P2M(&m, &c);
	assert(fabs(m.m - 4.5) < EPS && fabs(m.x - 14.75/4.5) < EPS && fabs(m.y - 19.25/4.5) < EPS);

	freeCell(&c);

	printf("OK\n\n");
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

	printf("OK\n\n");
}

void testDistance()
{
	printf("Testing cell-particle distance computation: ");

	Cell c;
	Multipole m;
	srand(42);
	initCell(&c, 10000, 10000, 1e30, 1e32, 0, 1e17, 0, 1e17);

	initMultipole(&m, 1.0, -1e17, 2e17, 0.0, 0.0, 0.0, 0.0);
	assert(fabs(distance(&m, &c) - sqrt(1e17*1e17 + 1e17*1e17)) < EPS);
	initMultipole(&m, 1.0, 0.5e17, 2e17, 0.0, 0.0, 0.0, 0.0);
	assert(fabs(distance(&m, &c) - 1e17) < EPS);
	initMultipole(&m, 1.0, 2e17, 2e17, 0.0, 0.0, 0.0, 0.0);
	assert(fabs(distance(&m, &c) - sqrt(1e17*1e17 + 1e17*1e17)) < EPS);

	initMultipole(&m, 1.0, -1e17, 0.5e17, 0.0, 0.0, 0.0, 0.0);
	assert(fabs(distance(&m, &c) - 1e17) < EPS);
	initMultipole(&m, 1.0, 0.5e17, 0.5e17, 0.0, 0.0, 0.0, 0.0);
	assert(fabs(distance(&m, &c) - 0) < EPS);
	initMultipole(&m, 1.0, 2e17, 0.5e17, 0.0, 0.0, 0.0, 0.0);
	assert(fabs(distance(&m, &c) - 1e17) < EPS);

	initMultipole(&m, 1.0, -1e17, -1e17, 0.0, 0.0, 0.0, 0.0);
	assert(fabs(distance(&m, &c) - sqrt(1e17*1e17 + 1e17*1e17)) < EPS);
	initMultipole(&m, 1.0, 0.5e17, -1e17, 0.0, 0.0, 0.0, 0.0);
	assert(fabs(distance(&m, &c) - 1e17) < EPS);
	initMultipole(&m, 1.0, 2e17, -1e17, 0.0, 0.0, 0.0, 0.0);
	assert(fabs(distance(&m, &c) - sqrt(1e17*1e17 + 1e17*1e17)) < EPS);
	
	printf("OK\n\n");
}

// TUNING TEST

void tuneM2P()
{
	printf("Tuning M2P. Force exerted by a cell on another cell of 10k particles. "
		   "M2P vs P2P_extRef.\n");

	
	Cell c0, c1, c2;
	Multipole m2;
	WorkingVecs wv;
	initWorkingVecs(&wv);
	
	double maxRelativeError, meanRelativeError, d;
	double l = 1e17;

	printf("Influence of the MAC, 10k particles per cell:\n"
		   "   l/d   | max relative error | average relative error:\n");
	for (double i = 1e17; i < 4e17; i += 0.25e17)
	{
		srand(42);
		initCell(&c0, 10000, 10000, 1.0, 100.0, 0, 1e17, 0, 1e17);
		srand(42);
		initCell(&c1, 10000, 10000, 1.0, 100.0, 0, 1e17, 0, 1e17);

		initCell(&c2, 10000, 10000, 1, 100.0, i, i+1e17, i, i+1e17);
		P2M(&m2, &c2);

		P2P_extRef(&c2, &c0);
		M2P(&m2, &c1, &wv);
		
		computeRelativeErrors(&c0, 1, &c1, &maxRelativeError, &meanRelativeError);
		d = distance(&m2, &c0);
		printf("%f |    %e    |       %e\n", l/d, maxRelativeError, meanRelativeError);
		
		freeCell(&c0);
		freeCell(&c1);
		freeCell(&c2);
	}

	printf("Influence of the nb of particles in the cells, l/d = 0.707:\n"
		   "nbParticles   | max relative error | average relative error:\n");
	for (int N = 100; N < 25000; N *=2)
	{
		srand(42);
		initCell(&c0, N, N, 1.0, 100.0, 0, 1e17, 0, 1e17);
		srand(42);
		initCell(&c1, N, N, 1.0, 100.0, 0, 1e17, 0, 1e17);

		initCell(&c2, N, N, 1, 100.0, 1.5e17, 2.5e17, 1.5e17, 2.5e17);
		P2M(&m2, &c2);
		
		P2P_extRef(&c2, &c0);
		M2P(&m2, &c1, &wv);
		computeRelativeErrors(&c0, 1, &c1, &maxRelativeError, &meanRelativeError);
		d = distance(&m2, &c0);
		printf("    %*d     |    %e    |       %e\n", 5, N, maxRelativeError, meanRelativeError);
		
		freeCell(&c0);
		freeCell(&c1);
		freeCell(&c2);
	}

	freeWorkingVecs(&wv);
}



// REGRESSION TESTS

void testPonP()
{
	printf("Regression test 1, pOnP vs pOnP_ref (10k particles):\n");

	Cell c1, c2;
	srand(42);
	initCell(&c1, 10000, 10000, 1e30, 1e32, 0, 1e17, 0, 1e17);
	srand(42);
	initCell(&c2, 10000, 10000, 1e30, 1e32, 0, 1e17, 0, 1e17);

	ponP_ref(1e31, 0.5e17, 0.5e17, c1.nbParticles, c1.m, c1.x, c1.y, c1.fx, c1.fy); 
	
	WorkingVecs wv;
	initWorkingVecs(&wv);
	ponP(1e31, 0.5e17, 0.5e17, c2.nbParticles, c2.m, c2.x, c2.y, c2.fx, c2.fy, &wv); 
	
	double maxRelativeError;
	double meanRelativeError;
	computeRelativeErrors(&c1, 1, &c2, &maxRelativeError, &meanRelativeError);

	printf("max relative error = %e, average relative error = %e\n\n", maxRelativeError, meanRelativeError);

	freeCell(&c1);
	freeCell(&c2);
	freeWorkingVecs(&wv);
}

void testP2Pin()
{
	printf("Regression test 2, P2P_in vs P2P_inRef (10k particles):\n");

	Cell c1, c2;
	srand(42);
	initCell(&c1, 10000, 10000, 1e30, 1e32, 0, 1e17, 0, 1e17);
	srand(42);
	initCell(&c2, 10000, 10000, 1e30, 1e32, 0, 1e17, 0, 1e17);

	P2P_inRef(&c1);

	WorkingVecs wv;
	initWorkingVecs(&wv);
	P2P_in(&c2, &wv);
	
	double maxRelativeError;
	double meanRelativeError;
	computeRelativeErrors(&c1, 1, &c2, &maxRelativeError, &meanRelativeError);

	printf("max relative error = %e, average relative error = %e\n\n", maxRelativeError, meanRelativeError);

	freeCell(&c1);
	freeCell(&c2);
	freeWorkingVecs(&wv);
}

void testP2Pext()
{
	printf("Regression test 3, P2P_ext vs P2P_extRef (10k particles):\n");

	Cell c1, c2, c3;
	srand(42);
	initCell(&c1, 10000, 10000, 1e30, 1e32, 0, 1e17, 0, 1e17);
	srand(42);
	initCell(&c2, 10000, 10000, 1e30, 1e32, 0, 1e17, 0, 1e17);
	
	initCell(&c3, 10000, 10000, 1e30, 1e32, 0, 1e17, 0, 1e17);

	P2P_extRef(&c3, &c1);

	WorkingVecs wv;
	initWorkingVecs(&wv);
	P2P_ext(&c3, &c2, &wv);
	
	double maxRelativeError;
	double meanRelativeError;
	computeRelativeErrors(&c1, 1, &c2, &maxRelativeError, &meanRelativeError);

	printf("max relative error = %e, average relative error = %e\n\n", maxRelativeError, meanRelativeError);

	freeCell(&c1);
	freeCell(&c2);
	freeWorkingVecs(&wv);
}

void testBarnesHut()
{
	srand(42);
	printf("Regression test 4, comparing results of a Barnes-Hut like method on a 4x4 grid" 
		   "(10 k particles) VS P2P_inRef on one global cell:\n");

	Cell cMerged;
	Quadtree qt;
	initQuadtree(&qt, 3, 10000, 10000, 1e30, 1e32, 0, 1e17, 0, 1e17);
	mergeCell(&cMerged, qt.nbCells, qt.cells);

	P2P_inRef(&cMerged);

	computeMultipoles(&qt);
	computeForces(&qt);	

	double maxRelativeError;
	double meanRelativeError;
	computeRelativeErrors(&cMerged, qt.nbCells, qt.cells, &maxRelativeError, &meanRelativeError);

	printf("max relative error = %e, average relative error = %e\n\n", maxRelativeError, meanRelativeError);

	freeCell(&cMerged);
	freeQuadtree(&qt);
}

void devBench()
{
	printf("Dev bench: ");
	Quadtree qt;
	Cell cMerged;

  	struct timeval start, stop;

	srand(42);
  	gettimeofday(&start, NULL);
	initQuadtree(&qt, 4, 200000, 200000, 1e30, 1e32, 0, 1e17, 0, 1e17);
  	gettimeofday(&stop, NULL);
	mergeCell(&cMerged, qt.nbCells, qt.cells);
	printf("Particles generation: %lf\n", stop.tv_sec - start.tv_sec + 0.000001 * (stop.tv_usec - start.tv_usec));

	gettimeofday(&start, NULL);
	computeMultipoles(&qt);
	computeForces(&qt);	
  	gettimeofday(&stop, NULL);
	printf("Interactions computation w/ barnes-hut like method: %lf\n", stop.tv_sec - start.tv_sec + 0.000001 * (stop.tv_usec - start.tv_usec));

	gettimeofday(&start, NULL);
	P2P_inRef(&cMerged);
  	gettimeofday(&stop, NULL);
	printf("Interactions computation w/ naive method: %lf\n", stop.tv_sec - start.tv_sec + 0.000001 * (stop.tv_usec - start.tv_usec));
	
	freeQuadtree(&qt);
	freeCell(&cMerged);
	printf("done\n");
}


int main(int argc, char const *argv[])
{
	// testP2M();
	// testInitQuadtree();
	// testDistance();

	// tuneM2P();
	
	// testPonP();
	// testP2Pin();
	// testP2Pext();
	// testBarnesHut();

	devBench();
	return EXIT_SUCCESS;
}