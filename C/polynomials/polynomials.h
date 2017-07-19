#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

#include "maths.h"

/*STRUCTURES*/

typedef struct Pol{
	int deg; /*degree of the polynomial*/
	double *coef; /*coefficients of the polynomials (from least to greatest)*/
} Pol;

/*PROTOTYPES*/

/*polynomials manipulation*/

void create_pol(Pol *pol, double deg, double *coef);
	/*Associates pol to the polynomial of degree def and coefficients coef.*/

void init_pol(Pol *pol, int len);
	/*Initialises the array of polynomials pol of length len to 0.*/

void free_mat_pol(Pol *mat, int len);
	/*Frees the polynomials in the matrix mat of length len.*/

void equal_pol(Pol *pol, Pol *pol_);
	/*Associates pol to pol_.*/

/*polynomials mathematical operations*/
/*the polynomials manipulated here can appear in the results and argument of the functions,
therefore all the polynomials passed to these functions have to be INITIALISED*/

void pol_prod(Pol *pol, Pol *pol1, Pol *pol2);
	/*Associates pol to the product of pol1 and pol2.*/

void pol_2prod(Pol *pol, Pol *pol1, Pol *pol2, Pol *pol3);
	/*Associates pol to the product of pol1, pol2 and pol3.*/

void pol_sum(Pol *pol, Pol *pol1, Pol *pol2);
	/*Associates pol to the sum of pol1 and pol2.*/

void pol_diff(Pol *pol, Pol *pol1, Pol *pol2);
	/*Associates pol to the difference of pol1 and pol2.*/

void scal_pol_prod(Pol *pol_, double scal, Pol *pol);
	/*Associates pol_ to scal * pol.*/

/*matrix of polynomials manipulation*/
/*the polynomials manipulated here can not appear in the results and argument of the functions,
therefore all the polynomials passed to these functions have to be UNINITIALISED*/

void scal2pol(Pol *pol, double *scal, int len);
	/*Associates the list of polynomials pol of length len to the list of scalars
	scal of length len.*/

void scal_mat_pol(Pol *mat, Pol *pol, int size);
	/*Associates mat to the scalar matrix of size size*size whose diagonal values
	are the polynomial pol.*/

void mat_pol_prod(Pol *pol, Pol *pol1, Pol *pol2, int i1, int j1, int j2);
	/*Associates to pol the matrix product of pol1 and pol2, with i1 the number of lines of pol1 and
	j1 and j2 the number of columns of pol1 and pol2.
	Notice that the number of lines of pol2 have to be i1.*/

void det2_mat_pol(Pol *pol, Pol *mat);
	/*Associates pol to the determinant of the 2*2 matrix of polynomials mat.*/

void det3_mat_pol(Pol *pol, Pol *mat);
	/*Associates pol to the determinant of the 3*3 matrix of polynomials mat.
	This function uses the Sarrus rule.*/

void adj_mat3_pol(Pol *adj, Pol *mat);
	/*Associates to adj the ajoint matrix of mat of size 3*3.*/

/*polynomial functions manipulation*/

double eval_pol(Pol *pol, double x);
	/*Evaluates the polynomial pol in x using the Horner's rule.*/

void der_pol(Pol *dpol, Pol *pol);
	/*Associates dpol to the derivative of pol.*/

double root_Newton_pol(Pol *pol, double start, double epsilon);
	/*Returns a root of the polynomial pol closed to start, found with the Newton's method with a
	precision epsilon, in the interval ]0,1[.*/

double root_Haley_pol(Pol *pol, double start, double epsilon);
	/*Returns a root of the polynomial pol closed to start, found with the Haley's method with a
	precision epsilon, in the interval ]0,1[.*/

#endif
