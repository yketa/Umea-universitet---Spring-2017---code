#ifndef ELLIPSOIDS_H
#define ELLIPSOIDS_H

#include "param.h"
#include "maths.h"

/*STRUCTURES*/

typedef struct Overlap{
	/*DESCRIBES THE INTERACTIONS BETWEEN PARTICLES*/
	/*defined per couple of particles*/
	/*Communication between functions*/
	int overlapping; /*boolean ; 1 = overlapping / 0 = not overlapping*/
	/*Variables for the computation of forces*/
	double mu; /*resclaing factor for the ellipsoids to be externally tangent*/
	/*refreshed only if overlap.overlapping*/
	double r_c[3]; /*contact point*/
	double r_c1[3]; /*contact point on particle 1*/
	double r_c2[3]; /*contact point on particle 2*/
	double v1c2[3]; /*local velocity of particle 1 at its point of contact with particle 2*/
	double v2c1[3]; /*local velocity of particle 2 at its point of contact with particle 1*/
	/*variables for measurements*/
	double r_c1_in1[3]; /*coordinates of the contact point on the surface of particle 1 in particle 1 frame*/
	double r_c2_in2[3]; /*coordinates of the contact point on the surface of particle 2 in particle 2 frame*/
	/*error checking*/
	double err; /*error on the determination of the root of polynomial H_{AB}*/
} Overlap;

typedef struct SizeEllipsoid{
	/*DESCRIBES THE PHYSICAL PARAMETERS OF THE PARTICLES*/
	/*chosen*/
	double m; /*mass of the ellipsoid*/
	double R; /*nominal radius*/
	double major; /*coefficient of the major axis*/
	double minor; /*coefficient of the minor axis*/
	/*inferred*/
	double ax[3]; /*semi-axes of the ellipsoid*/
	double I[3]; /*moments of inertia of the ellipsoid*/
	int axMax; /*index of the longest axis*/
} SizeEllipsoid;

typedef struct Ellipsoid{
	/*DESCRIBES THE PARTICLES AND THEIR MOVEMENTS*/
	/*Variables to measure*/
	double r[3]; /*centre of the ellipsoid*/
	double v[3]; /*velocity of the elliposoid*/
	double q[4]; /*quaternion of the ellipsoid*/
	double w[3]; /*rotation vector of the particle in the body frame*/
	/*Description of the ellipsoid*/
	int size; /*index of the size of the ellipsoid*/
	double Q[9]; /*rotation matrix associated to q*/
} Ellipsoid;

typedef struct Forces{
	/*DESCRIBES THE EVOLUTION OF THE MOVEMENTS OF THE PARTICLES*/
	/*Forces*/
	double el[3]; /*elastic forces*/
	double dis[3]; /*dissipative forces*/
	/*Moments*/
	double M[3]; /*moments of forces*/
	/*Quaternions*/
	double wdot[3]; /*first derivative of the rotation vector in the body frame*/
} Forces;

/*PROTOTYPES*/

void rotation_matrix(Ellipsoid *ellipsoid);
	/*Updates the rotation matrix associated to ellipsoid->q.*/

double belonging_function(Ellipsoid *ellipsoid, SizeEllipsoid *sizes, double *vec);
	/*Returns the belonging fonction of ellipsoid, with semi-axes ax, evaluated in vec.*/

void omega(double *w, Ellipsoid *ellipsoid);
	/*Associates w to the rotation vector of ellipsoid.*/

void qdot(double *qp, Ellipsoid *ellipsoid);
	/*Associates qp to the first derivative of the quaternion of ellipsoid.*/

void wdot(Ellipsoid *ellipsoid, SizeEllipsoid *sizes, Forces *forces);
	/*Associates to forces->wdot the first derivative of the rotation vector of the ellipsoid.*/

void n_surface_vec(double *n, Ellipsoid *ellipsoid, SizeEllipsoid *sizes, double *vec);
	/*Associates n to the non-unit surface vector of ellipsoid, whose semi-axes are ax,
	at vec.*/

void inv_red_bel(double *Binv, Ellipsoid *ellipsoids, SizeEllipsoid *sizes);
	/*Associates Binv to the inverse of the reduced belonging matrix of ellipsoid.*/

void adj_y_ab(double *adj, double *Bainv, double *Bbinv);
	/*Associates to adj the adjoint of the matrix Y_{AB} defined by the inverses of
	the reduced belonging matrix of the ellipsoids, Bainv and Bbinv.
	Coefficients calculated with Mathematica.*/

void p_ab(double *pab, double *adj, double *rab);
	/*Associates to pab the coefficients of the polynomial P_{AB} corresponding to
	the ajoint adj of Y_{AB} and ellipsoids whose centers are separated by vector
	rab = \vec{r_B} - \vec{r_A}.*/

void q_ab(double *det, double *Bainv, double *Bbinv);
	/*Associated to det the coefficients of the determinant of the matrix Y_{AB}
	defined the inverses of the reduced belonging matrix of the ellipsoids, Bainv
	and Bbinv.
	Coefficients calculated with Mathematica.*/

double root_Haley_pol_hAB(double *hab, double *dhab, double start, double epsilon);
	/*Returns the root of the polynomial H_{AB} closed to start, found with the Haley's method with a
	precision epsilon, in the interval ]0,1[.*/

double eval_pol(double *coef, int deg, double x);
	/*Evaluates the polynomial of coefficients coef, degree deg
	in x using the Horner's rule.*/

void contact(Par *par, Overlap *overlap, SizeEllipsoid *sizes,
	Ellipsoid *ellipsoid1, Ellipsoid *ellipsoid2);
	/*Updates the contact points between ellipsoid1 and ellipsoid2.
	Algorithm described in subsection 2.5.1 — based on Perram and Wertheim, J. of Comp. Phys., 1995*/

#endif
