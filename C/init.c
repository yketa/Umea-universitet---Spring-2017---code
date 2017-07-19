/*REFERENCE: https://fr.sharelatex.com/project/59156f86e054ce7b2148b2b2*/

/* This file contains the number of ellipsoids, their sizes, and the initial conditions
for the simulation run by main.c.*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "init.h"

// ---------- EDIT ----------


// PHYSICAL PARAMETERS

double ke = 1; // constant of elasticity
double kd = 1; // constant of dissipation*/

/*SIZES*/

int nSizes = 1; // number of sizes of ellipsoids

double m[] = {1}; // masses — as a 1xnSizes matrix
double R[] = {0.5}; // nominal radii — as a 1xnSizes matrix
double major[] = {1.1}; // coefficients of the major axis — as a 1xnSizes matrix
double minor[] = {0.9534625892455922}; // coefficients of the minor axis — as a 1xnSizes matrix

// ELLIPSOIDS

int nEllipsoids = 4; // number of ellipsoids

int eSizes[] = {0,0,0,0}; // index of the sizes of the ellipsoids

double r[] = {0,0,0,3,0.5,0,-3,0.5,0,0,0,4}; // initial positions — as a 3xnEllipsoids matrix
double v[] = {0,0,0,-0.01,0,0,0.02,0,0,0,0,-0.03}; // initial velocities — as a 3xnEllipsoids matrix
double q[] = {0,0,0,1,0,0,0.7071067811865476,0.7071067811865476,0,0.7071067811865476,0,0.7071067811865476,0,0,0,1}; // initial quaternions — as a 4xnEllipsoids matrix
double w[] = {0,0,0,0,0,0,0,0,0,0,0,0}; // initial angular velocities — as a 3xnEllipsoids matrix

// INTEGRATION PARAMETERS

double dt = 0.05; // time step
int Niter = 10000; // numer of iterations

double epsilon_dist = 1e-7; // precision on distances
double epsilon_ang = 1e-3; // precision on angles
double epsilon_mu = 0.5; // half width of the interval of search of the rescaling factor
double epsilon = 1e-12; // precision


// ---------- DO NOT EDIT ----------


// INITIALISATION

void init_struct(Par *par, SizeEllipsoid *sizes, Ellipsoid *ellipsoids){
	// Initialisation of the paramaters, sizes and ellipsoids structures according
	// to the parameters in init.c.
	
	// PHYSICAL PARAMETERS

	double axinv2[3];

	for (int s = 0; s < nSizes; s++){
		(sizes + s)->m = m[s];

		// semi-axes
		(sizes + s)->R = R[s];
		(sizes + s)->major = major[s];
		(sizes + s)->minor = minor[s];

		(sizes + s)->ax[0] = (sizes + s)->major*(sizes + s)->R;
		(sizes + s)->ax[1] = (sizes + s)->minor*(sizes + s)->R;
		(sizes + s)->ax[2] = (sizes + s)->minor*(sizes + s)->R;

		for (int coord = 0; coord < 3; coord++){
			axinv2[coord] = 1/pow((sizes + s)->ax[coord],2);
		}
		// (sizes + s)->g = 1/axinv2[max(axinv2,3)]; // gamma factor of ever-belonging balls

		(sizes + s)->axMax = max((sizes + s)->ax,3); // index of the longest axis [OBSOLETE]

		for (int coord = 0; coord < 3; coord++){
			(sizes + s)->I[coord] = 0.2*(sizes + s)->m*(pow((sizes + s)->ax[(coord + 1)%3],2) + pow((sizes + s)->ax[(coord + 2)%3],2)); // moments of inertia
		}
	}

	par->ke = ke; // constant of elasticity
	par->kd = kd; // constant of dissipation

	// INITIAL PARAMETERS

	for (int el = 0; el < nEllipsoids; el++){
		for (int coord = 0; coord < 3; coord++){
			(ellipsoids + el)->r[coord] = r[el*3 + coord]; // initial positions
			(ellipsoids + el)->v[coord] = v[el*3 + coord]; // initial velocities
			(ellipsoids + el)->q[coord] = q[el*4 + coord]; // initial quaternions
			(ellipsoids + el)->w[coord] = w[el*3 + coord]; // initial first derivatives of quaternions
		}
		(ellipsoids + el)->q[3] = q[el*4 + 3];

		(ellipsoids + el)->size = eSizes[el]; // index of the size of the ellipsoid
		rotation_matrix(ellipsoids + el); // rotation matrix
	}

	// INTEGRATION PARAMETERS

	par->epsilon_dist2 = pow(epsilon_dist,2); // precision on distances
	par->epsilon_ang = epsilon_ang; // precision on angles
	par->epsilon_mu = epsilon_mu; // half width of the interval of search of the rescaling factor
	par->epsilon = epsilon; // precision

	par->Niter = Niter; // number of iterations

}
