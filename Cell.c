#include "Cell.h"
#include "Particle.h"

#include <malloc.h>
#include <stdlib.h>

#define NB_SUBCENTERS 4


// Initialize a cell with a random number of particles between
// nbParticlesMin and nbParticlesMax, of random mass between mMin and mMax, and random 
// positions (x,y) with x between xMin and xMax and y between yMin and yMax
void initCell(Cell *c, int nbParticlesMin, int nbParticlesMax, double mMin, double mMax, double xMin, double xMax, double yMin, double yMax)
{
	c->nbParticles = (rand()/(float)(RAND_MAX)) * (nbParticlesMax - nbParticlesMin) + nbParticlesMin;
	c->particles = (Particle *) malloc(c->nbParticles * sizeof(Particle));

	int m, x, y;
	for (int i = 0; i < c->nbParticles; i++)
	{
		m = (rand()/(float)(RAND_MAX)) * (mMax - mMin) + mMin;
		x = (rand()/(float)(RAND_MAX)) * (xMax - xMin) + xMin;
		y = (rand()/(float)(RAND_MAX)) * (yMax - yMin) + yMin;
		initParticle(c->particles + i, m, x, y);
	}
}

// Apply the forces of the particles in c on each other
void P2P_in(Cell *c)
{

	for (int i = 0; i < c->nbParticles; i++)
		for (int j = 0; j < c->nbParticles; j++)
			if (i != j)
				PonP(c->particles + i, c->particles + j);				
}

// Apply the forces of the particles in c1 on the particles in c2
void P2P_ext(Cell *c1, Cell *c2)
{
	for (int i = 0; i < c1->nbParticles; i++)
		for (int j = 0; j < c2->nbParticles; j++)
				PonP(c1->particles + i, c2->particles + j);				
}

// Initialize m as center of mass of particles in cell c
void P2M(Particle *m, Cell *c)
{
	m->m = m->x = m->x = m->fx = m->fy = 0;

	for (int i = 0; i < c->nbParticles; i++)
	{
		m->m += c->particles[i].m;
		m->x += c->particles[i].m * c->particles[i].x;
		m->y += c->particles[i].m * c->particles[i].y;
	}

	m->x /= m->m;
	m->y /= m->m;
}

// Apply the forces of a center of mass m on particles in cell c
void M2P(Particle *m, Cell *c)
{
	for (int i = 0; i < c->nbParticles; i++)
		PonP(m, c->particles + i);				
}

// Initialize m as center of mass of centers of mass subM[0], subM[1], subM[2] subM[3]
void M2M(Particle *m, Particle *subM)
{
	m->m = m->x = m->x = m->fx = m->fy = 0;

	for (int i = 0; i < NB_SUBCENTERS; i++)
	{
		m->m += subM[i].m;
		m->x += subM[i].m * subM[i].x;
		m->y += subM[i].m * subM[i].y;
	}

	m->x /= m->m;
	m->y /= m->m;	
}


// Release ressources associated with the cell c
void freeCell(Cell *c)
{
	free(c->particles);
}