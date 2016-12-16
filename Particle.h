#ifndef PARTICLE_H
#define PARTICLE_H

#include <malloc.h>

// #define G 6.67e-11 // m^3.kg^-1.s^-2

struct Particle
{
  double m; // mass
  double x; // x coordinate in the plane
  double y; // y coordinate in the plane

  double vx; // speed along x-axis
  double vy; // speed along y-axis

  double fx; // gravitationnal force exerted on the particle
  double fy; // gravitationnal force exerted on the particle

} typedef Particle;

Particle *initParticle(double m, double x, double y, double vx, double vy);

void freeParticle(Particle *p);



#endif
