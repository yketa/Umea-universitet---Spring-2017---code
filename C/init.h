#ifndef INIT_H
#define INIT_H

#include "ellipsoids.h"
#include "maths.h"

/*PHYSICAL PARAMETERS*/

extern double ke; /*constant of elasticity*/
extern double kd; /*constant of dissipation*/

/*SIZES*/

extern int nSizes; /*number of sizes of ellipsoids*/

extern double m[]; /*masses — as a 1xnSizes matrix*/
extern double R[]; /*nominal radii — as a 1xnSizes matrix*/
extern double major[]; /*coefficients of the major axis — as a 1xnSizes matrix*/
extern double minor[]; /*coefficients of the minor axis — as a 1xnSizes matrix*/

/*ELLIPSOIDS*/

extern int nEllipsoids; /*number of ellipsoids*/

extern int eSizes[]; /*index of the sizes of the ellipsoids*/

extern double r[]; /*initial positions — as a 3xnEllipsoids matrix*/
extern double v[]; /*initial velocities — as a 3xnEllipsoids matrix*/
extern double q[]; /*initial quaternions — as a 4xnEllipsoids matrix*/
extern double w[]; /*initial angular velocities — as a 3xnEllipsoids matrix*/

/*INTEGRATION PARAMETERS*/

extern double dt; /*time step*/
extern int Niter; /*numer of iterations*/

extern double epsilon_dist; /*precision on distances*/
extern double epsilon_ang; /*precision on angles*/
extern double epsilon_mu; /*half width of the interval of search of the rescaling factor*/
extern double epsilon; /*precision*/

/*PROTOTYPES*/

void init_struct(Par *par, SizeEllipsoid *sizes, Ellipsoid *ellipsoids);
	/*Initialisation of the paramaters, sizes and ellipsoids structures according
	to the parameters in init.c.*/

#endif
