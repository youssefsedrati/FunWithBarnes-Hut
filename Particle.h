#ifndef PARTICLE_H
#define PARTICLE_H

#include <malloc.h>

// Real particle or virtual particle (multipole, fx & fy are then not used)
struct Particle
{
  double m; // mass
  double x; // x coordinate in the plane
  double y; // y coordinate in the plane

  double fx; // force along x-axis
  double fy; // force along y-axis

} 
typedef Particle;

// Init a particle of mass m, and position (x,y)
void initParticle(Particle *p, double m, double x, double y);

// Apply the force of p1 on p2
void PonP(Particle *p1, Particle *p2);


#endif
