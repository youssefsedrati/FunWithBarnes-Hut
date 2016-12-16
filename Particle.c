#include "Particle.h"

Particle *initParticle(double m, double x, double y, double vx, double vy)
{
  Particle *p = (Particle *) malloc(sizeof(Particle));
  p->m = m;
  p->x = x;
  p->y = y;
  p->vx = vx;
  p->vy = vy;
  return p;
}

void freeParticle(Particle *p)
{
  free(p);
}
