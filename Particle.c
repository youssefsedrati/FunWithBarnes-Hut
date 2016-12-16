#include "Particle.h"

#include <math.h>

#define G 6.67e-11 // m^3.kg^-1.s^-2

// Init a particle of mass m, and position (x,y)
void initParticle(Particle *p, double m, double x, double y)
{
  p->m = m;
  p->x = x;
  p->y = y;
  p->fx = 0;
  p->fy = 0;
}

// Apply the force of p1 on p2
void PonP(Particle *p1, Particle *p2)
{
  double dx = p2->x - p1->x;
  double dy = p2->y - p1->y;
  double d = sqrt(dx*dx + dy*dy);
  double r = G * p1->m * p2->m / (d*d*d);
  p2->fx -= r * dx;
  p2->fy -= r * dy;
}

