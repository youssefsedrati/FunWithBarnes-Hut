#include "Cell.h"
#include "Particle.h"

#include <malloc.h>
#include <stdlib.h>
#include <string.h>

#define NB_SUBCENTERS 4


// Initialize a cell with a random number of particles between
// nbParticlesMin and nbParticlesMax, of random mass between mMin and mMax, and random 
// positions (x,y) with x between xMin and xMax and y between yMin and yMax
void initCell(Cell *c, int nbParticlesMin, int nbParticlesMax, double mMin, double mMax, double xMin, double xMax, double yMin, double yMax)
{
	c->xMin = xMin;
	c->xMax = xMax;
	c->yMin = yMin;
	c->yMax = yMax;
	c->nbParticles = (rand()/(float)(RAND_MAX)) * (nbParticlesMax - nbParticlesMin) + nbParticlesMin;
	c->particles = (Particle *) malloc(c->nbParticles * sizeof(Particle));

	double m, x, y;
	for (int i = 0; i < c->nbParticles; i++)
	{
		m = (rand()/(float)(RAND_MAX)) * (mMax - mMin) + mMin;
		x = (rand()/(float)(RAND_MAX)) * (xMax - xMin) + xMin;
		y = (rand()/(float)(RAND_MAX)) * (yMax - yMin) + yMin;
		initParticle(c->particles + i, m, x, y);
	}
}

// Merge cells in cMerged
void mergeCell(Cell *cMerged, int nbCells, Cell *cells)
{
	cMerged->nbParticles = 0;
	for (int i = 0; i < nbCells; i++)
		cMerged->nbParticles += cells[i].nbParticles;

	cMerged->particles = (Particle *) malloc(cMerged->nbParticles * sizeof(Particle));

	int p = 0;
	for (int i = 0; i < nbCells; i++)
	{
		memcpy(cMerged->particles + p, cells[i].particles, cells[i].nbParticles * sizeof(Particle));
		p += cells[i].nbParticles;
	}
}

// Apply the forces of the particles in c on each other
void P2P_in(Cell *c)
{
	for (int i = 0; i < c->nbParticles; i++)
		for (int j = 0; j < i; j++)
		{
			PonP(c->particles + i, c->particles + j);
			PonP(c->particles + j, c->particles + i);				
		}
}

// Apply the forces of the particles in c1 on the particles in c2
void P2P_ext(Cell *c1, Cell *c2)
{
	for (int i = 0; i < c1->nbParticles; i++)
		for (int j = 0; j < c2->nbParticles; j++)
				PonP(c1->particles + i, c2->particles + j);				
}

// Initialize cm as center of mass of particles in cell c
void P2M(Particle *cm, Cell *c)
{
	cm->m = cm->x = cm->y = cm->fx = cm->fy = 0;

	for (int i = 0; i < c->nbParticles; i++)
	{
		cm->m += c->particles[i].m;
		cm->x += c->particles[i].m * c->particles[i].x;
		cm->y += c->particles[i].m * c->particles[i].y;
	}

	cm->x /= cm->m;
	cm->y /= cm->m;
}

// Apply the forces of a center of mass cm on particles in cell c
void M2P(Particle *cm, Cell *c)
{
	for (int i = 0; i < c->nbParticles; i++)
		PonP(cm, c->particles + i);				
}

// Initialize cm as center of mass of centers of mass subCms[0], subCms[1], subCms[2] subCms[3]
void M2M(Particle *cm, Particle *subCms)
{
	cm->m = cm->x = cm->y = cm->fx = cm->fy = 0;

	for (int i = 0; i < NB_SUBCENTERS; i++)
	{
		cm->m += subCms[i].m;
		cm->x += subCms[i].m * subCms[i].x;
		cm->y += subCms[i].m * subCms[i].y;
	}

	cm->x /= cm->m;
	cm->y /= cm->m;	
}


// Release ressources associated with the cell c
void freeCell(Cell *c)
{
	free(c->particles);
}