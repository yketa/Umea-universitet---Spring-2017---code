#ifndef MEASURE_H
#define MEASURE_H

#include "init.h"
#include "maths.h"
#include "ellipsoids.h"

/*PROTOTYPES*/

FILE* init_measure(int argc, char *argv[], Par *par, SizeEllipsoid *sizes, Ellipsoid *ellipsoids, Forces *forces,
	Overlap *overlaps, double *p, double *L);
	/*Initialises the meausrement file.*/

void measure(FILE* file, Par *par, SizeEllipsoid *sizes, Ellipsoid *ellipsoids, Forces *forces, Overlap *overlaps,
	double *p, double *L);
	/*Prints in .csv format the time and the positions, velocities and angular velocities of the ellipsoids.*/

void momentum(Ellipsoid *ellipsoids, SizeEllipsoid *sizes, double *p, double *L);
	/*Associates to p and L the total translational and angular momentum on each axis
	on two consecutives time steps.*/

#endif
