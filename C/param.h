#ifndef PARAM_H
#define PARAM_H

/*STRUCTURES*/

typedef struct Par{
	/*Forces parameters*/
	double ke; /*constant of elasticity*/
	double kd; /*constant of dissipation*/
	/*Algorithms parameters*/
	double epsilon_dist2; /*precision on distances squared*/
	double epsilon_ang; /*precision on angles*/
	double epsilon_mu; /*half width of the interval of search of the rescaling factor*/
	double epsilon; /*precision*/
	/*Integration parameters*/
	int Niter; /*numer of iterations*/
	int iter; /*index of the current interation*/
} Par;

#endif
