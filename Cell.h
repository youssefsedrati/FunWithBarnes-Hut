#ifndef CELL_H
#define CELL_H

#include "Particle.h"

struct Cell
{
	double xMin;
	double xMax;
	double yMin;
	double yMax;
	int nbParticles;
	Particle *particles;
} 
typedef Cell;

// Initialize a cell with a random number of particles between
// nbParticlesMin and nbParticlesMax, of random mass between mMin and mMax, and random 
// positions (x,y) with x between xMin and xMax and y between yMin and yMax
void initCell(Cell *c, int nbParticlesMin, int nbParticlesMax, double mMin, double mMax, double xMin, double xMax, double yMin, double yMax);

// Merge cells in cMerged
void mergeCell(Cell *cMerged, int nbCells, Cell *cells);

// Apply the forces of the particles in c on each other
void P2P_in(Cell *c);

// Apply the forces of the particles in c1 on the particles in c2
void P2P_ext(Cell *c1, Cell *c2);

// Initialize m as center of mass of particles in cell c
void P2M(Particle *cm, Cell *c);

// Apply the forces of a center of mass m on particles in cell c
void M2P(Particle *cm, Cell *c);

// Initialize m as center of mass of centers of mass subM[0], subM[1], subM[2] subM[3]
void M2M(Particle *cm, Particle *subCms);

// Release ressources associated with the cell c
void freeCell(Cell *c);

#endif