#include "Multipole.h"

#include <math.h>

#define G 6.67e-11 // m^3.kg^-1.s^-2

// Initialize a multipole
void initMultipole(Multipole *mp, double m, double x, double y, 
				   double xMin, double xMax, double yMin, double yMax)

{
  mp->m = m;
  mp->x = x;
  mp->y = y;
  mp->xMin = xMin;
  mp->xMax = xMax;
  mp->yMin = yMin;
  mp->yMax = yMax;
}