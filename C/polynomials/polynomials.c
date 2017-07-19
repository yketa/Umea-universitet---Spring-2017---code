#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "polynomials.h"

/*polynomials manipulation*/

void create_pol(Pol *pol, double deg, double *coef){
	// Associates pol to the polynomial of degree def and coefficients coef.

	pol->deg = deg; // degree of the polynomial

	pol->coef = malloc(sizeof(double) * (pol->deg + 1));
	for (int co = 0; co <= pol->deg; co++){
		pol->coef[co] = coef[co]; // cofficients of the polynomial
	}
}

void init_pol(Pol *pol, int len){
	// Initialises the array of polynomials pol of length len to 0.

	for (int p = 0; p < len; p++){
		(pol + p)->deg = 0;

		(pol + p)->coef = malloc(sizeof(double));
		(pol + p)->coef[0] = 0;
	}
}

void free_mat_pol(Pol *mat, int len){
	// Frees the polynomials in the matrix mat of length len.

	for (int p = 0; p < len; p++){
		free((mat + p)->coef);
	}
}

void equal_pol(Pol *pol, Pol *pol_){
	// Associates pol to pol_.

	pol->deg = pol_->deg;

	pol->coef = malloc(sizeof(double) * (pol->deg + 1));
	equal(pol->coef,pol_->coef,pol->deg + 1);
}

/*polynomials mathematical operations*/
/*the polynomials manipulated here can appear in the results and argument of the functions,
therefore all the polynomials passed to these functions have to be INITIALISED*/

void pol_prod(Pol *pol, Pol *pol1, Pol *pol2){
	// Associates pol to the product of pol1 and pol2.

	double prod[pol1->deg + pol2->deg + 1];
	init(prod,pol1->deg + pol2->deg + 1); // initialisation of the coefficients
	for (int co1 = 0; co1 <= pol1->deg; co1++){
		for (int co2 = 0; co2 <= pol2->deg; co2++){
			prod[co1 + co2] += pol1->coef[co1] * pol2->coef[co2]; // coefficients of the polynomial
		}
	}

	pol->deg = pol1->deg + pol2->deg; // degree of the product
	free(pol->coef);
	pol->coef = malloc(sizeof(double) * (pol->deg + 1));
	equal(pol->coef,&prod[0],pol->deg + 1);
}

void pol_2prod(Pol *pol, Pol *pol1, Pol *pol2, Pol *pol3){
	// Associates pol to the product of pol1, pol2 and pol3.

	double prod[pol1->deg + pol2->deg + pol3->deg + 1];
	init(prod,pol1->deg + pol2->deg + pol3->deg + 1); // initialisation of the coefficients
	for (int co1 = 0; co1 <= pol1->deg; co1++){
		for (int co2 = 0; co2 <= pol2->deg; co2++){
			for (int co3 = 0; co3 <= pol3->deg; co3++){
				prod[co1 + co2 + co3] += pol1->coef[co1] * pol2->coef[co2] * pol3->coef[co3]; // coefficients of the polynomial
			}
		}
	}

	pol->deg = pol1->deg + pol2->deg + pol3->deg; // degree of the product
	free(pol->coef);
	pol->coef = malloc(sizeof(double) * (pol->deg + 1));
	equal(pol->coef,&prod[0],pol->deg + 1);
}

void pol_sum(Pol *pol, Pol *pol1, Pol *pol2){
	// Associates pol to the sum of pol1 and pol2.

	if (pol1->deg > pol2->deg){
		double sum[pol1->deg + 1];
		for (int co = 0; co <= pol2->deg; co++){
			sum[co] = pol1->coef[co] + pol2->coef[co]; // coefficients of the polynomial
		}
		for (int co = pol2->deg + 1; co <= pol1->deg; co++){
			sum[co] = pol1->coef[co]; // coefficients of the polynomial
		}

		pol->deg = pol1->deg; // degree of the polynomial
		free(pol->coef);
		pol->coef = malloc(sizeof(double) * (pol->deg + 1));
		equal(pol->coef,&sum[0],pol->deg + 1);
	}
	else {
		double sum[pol2->deg + 1];
		for (int co = 0; co <= pol1->deg; co++){
			sum[co] = pol1->coef[co] + pol2->coef[co]; // coefficients of the polynomial
		}
		for (int co = pol1->deg + 1; co <= pol2->deg; co++){
			sum[co] = pol2->coef[co]; // coefficients of the polynomial
		}

		pol->deg = pol2->deg; // degree of the polynomial
		free(pol->coef);
		pol->coef = malloc(sizeof(double) * (pol->deg + 1));
		equal(pol->coef,&sum[0],pol->deg + 1);
	}
}

void pol_diff(Pol *pol, Pol *pol1, Pol *pol2){
	// Associates pol to the difference of pol1 and pol2.

	if (pol1->deg > pol2->deg){
		double diff[pol1->deg + 1];
		for (int co = 0; co <= pol2->deg; co++){
			diff[co] = pol1->coef[co] - pol2->coef[co]; // coefficients of the polynomial
		}
		for (int co = pol2->deg + 1; co <= pol1->deg; co++){
			diff[co] = pol1->coef[co]; // coefficients of the polynomial
		}

		pol->deg = pol1->deg; // degree of the polynomial
		free(pol->coef);
		pol->coef = malloc(sizeof(double) * (pol->deg + 1));
		equal(pol->coef,&diff[0],pol->deg + 1);
	}
	else {
		double diff[pol2->deg + 1];
		for (int co = 0; co <= pol1->deg; co++){
			diff[co] = pol1->coef[co] - pol2->coef[co]; // coefficients of the polynomial
		}
		for (int co = pol1->deg + 1; co <= pol2->deg; co++){
			diff[co] = -pol2->coef[co]; // coefficients of the polynomial
		}

		pol->deg = pol2->deg; // degree of the polynomial
		free(pol->coef);
		pol->coef = malloc(sizeof(double) * (pol->deg + 1));
		equal(pol->coef,&diff[0],pol->deg + 1);
	}
}

void scal_pol_prod(Pol *pol_, double scal, Pol *pol){
	// Associates pol_ to scal * pol.

	double coef[pol->deg + 1]; // coefficients of the polynomial pol_
	for (int co = 0; co <= pol->deg; co++){
		coef[co] = scal * pol->coef[co]; // scalar multiplication
	}

	pol_->deg = pol->deg; // degree of the polynomial pol_
	free(pol_->coef);
	pol_->coef = malloc(sizeof(double) * (pol_->deg + 1));
	equal(pol_->coef,&coef[0],pol_->deg + 1);
}

/*matrix of polynomials manipulation*/
/*the polynomials manipulated here can not appear in the results and argument of the functions,
therefore all the polynomials passed to these functions have to be UNINITIALISED*/

void scal2pol(Pol *pol, double *scal, int len){
	// Associates the list of polynomials pol of length len to the list of scalars
	// scal of length len.

	for (int s = 0; s < len; s++){
		(pol + s)->deg = 0; // a scalar is a polynomial of degree 0

		(pol + s)->coef = malloc(sizeof(double)); // pol is uninitialised
		(pol + s)->coef[0] = *(scal + s); // scalar to polynomial
	}
}

void scal_mat_pol(Pol *mat, Pol *pol, int size){
	// Associates mat to the scalar matrix of size size*size whose diagonal values
	// are the polynomial pol.

	double coef[pol->deg + 1]; // coefficients of the polynomial pol
	equal(&coef[0],pol->coef,pol->deg + 1);

	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			if (i == j){
				// diagonal terms
				(mat + i*(size + 1))->deg = pol->deg;

				(mat + i*(size + 1))->coef = malloc(sizeof(double) * ((mat + i*(size + 1))->deg + 1)); // mat is uninitialised
				equal((mat + i*(size + 1))->coef,&coef[0],(mat + i*(size + 1))->deg + 1);
			}
			else {
				// non-diagonal terms
				init_pol(mat + (i * size + j),1); // mat is uninitialised
			}
		}
	}
}

void mat_pol_prod(Pol *pol, Pol *pol1, Pol *pol2, int i1, int j1, int j2){
	// Associates to pol the matrix product of pol1 and pol2, with i1 the number of lines of pol1 and
	// j1 and j2 the number of columns of pol1 and pol2.
	// Notice that the number of lines of pol2 have to be i1.

	Pol prod; // product of matrix coefficients
	init_pol(&prod,1);

	init_pol(pol,i1*j2); // pol is uninitialised
	for (int i = 0; i < i1; i++){
		for (int j = 0; j < j2; j++){
			for (int k = 0; k < j1; k++){
				pol_prod(&prod, pol1 + (i * j1 + k), pol2 + (k * j2 + j));
				pol_sum(pol + (i * j2 + j), pol + (i * j2 + j), &prod);
			}
		}
	}

	free(prod.coef);
}

void det2_mat_pol(Pol *pol, Pol *mat){
	// Associates pol to the determinant of the 2*2 matrix of polynomials mat.

	init_pol(pol,1); // pol is non-initialised

	pol_prod(pol,mat + 0,mat + 3); // first product in the determination of the determinant

	Pol prod_; // second product in the determination of the determinant
	init_pol(&prod_,1);
	pol_prod(&prod_,mat + 1,mat + 2);

	pol_diff(pol,pol,&prod_); // determinant

	free(prod_.coef);
}

void det3_mat_pol(Pol *pol, Pol *mat){
	// Associates pol to the determinant of the 3*3 matrix of polynomials mat.
	// This function uses the Sarrus rule.

	init_pol(pol,1); // pol is uninitialised

	Pol prod_; // second to sixth product in the determination of the determinant
	init_pol(&prod_,1);

	// first product
	pol_2prod(pol, mat + 0, mat + 4, mat + 8);

	// second product
	pol_2prod(&prod_, mat + 3, mat + 7, mat + 2);
	pol_sum(pol,pol,&prod_);

	// third product
	pol_2prod(&prod_, mat + 6, mat + 1, mat + 5);
	pol_sum(pol,pol,&prod_);

	// fourth product
	pol_2prod(&prod_, mat + 6, mat + 4, mat + 2);
	pol_diff(pol,pol,&prod_);

	// fifth product
	pol_2prod(&prod_, mat + 0, mat + 7, mat + 5);
	pol_diff(pol,pol,&prod_);

	// sixth product
	pol_2prod(&prod_, mat + 3, mat + 1, mat + 8);
	pol_diff(pol,pol,&prod_);

	free(prod_.coef);
}

void adj_mat3_pol(Pol *adj, Pol *mat){
	// Associates to adj the ajoint matrix of mat of size 3*3.

	Pol det2[4]; // matrix for which we need the determinant
	int count; // count of items in the minor matrix

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){

			// determination of the minor matrix
			count = 0;
			for (int i_ = 0; i_ < 3; i_++){
				for (int j_ = 0; j_ < 3; j_++){
					if (i_ != i && j_ != j){
						equal_pol(&det2[count],mat + (i_ * 3 + j_));
						count++;
					}
				}
			}

			// determination of the determinant
			det2_mat_pol(adj + (j*3 + i),&det2[0]);
			if ((i + j)%2 == 1){
				scal_pol_prod(adj + (j*3 + i),-1,adj + (j*3 + i));
			}

			// free det2
			free_mat_pol(&det2[0],4);

		}
	}
}

/*polynomial functions manipulation*/

double eval_pol(Pol *pol, double x){
	// Evaluates the polynomial pol in x using the Horner's rule.

	double result = pol->coef[pol->deg]; // a_n
	for (int c = pol->deg - 1; c >= 0; c--){
		result = pol->coef[c] + result * x; // a_{n-1} + b_n x
	}

	return result;
}

void der_pol(Pol *dpol, Pol *pol){
	// Associates dpol to the derivative of pol.

	if (pol->deg == 0){
		// pol is a constant polynomial
		dpol->deg = 0;

		free(dpol->coef);
		dpol->coef = malloc(sizeof(double));
		dpol->coef[0] = 0;
	}
	else{
		// pol is a non-constant polynomial

		double coef[pol->deg];
		for (int c = 0; c <= pol->deg - 1; c++){
			coef[c] = (c + 1)*pol->coef[c + 1]; // coefficients of the derivative of pol

		}

		dpol->deg = pol->deg - 1;
		free(dpol->coef);
		dpol->coef = malloc(sizeof(double) * (dpol->deg + 1));
		equal(dpol->coef,&coef[0],dpol->deg + 1);
	}
}

double root_Newton_pol(Pol *pol, double start, double epsilon){
	// Returns a root of the polynomial pol closed to start, found with the Newton's method with a
	// precision epsilon, in the interval ]0,1[.

	Pol dpol; // derivative of pol
	init_pol(&dpol,1);
	der_pol(&dpol,pol);

	while (pow(eval_pol(pol,start),2) > pow(epsilon,2) || start <= 0 || start >= 1){
		start -= eval_pol(pol,start)/eval_pol(&dpol,start);
	}

	free(dpol.coef);

	return start;
}

double root_Haley_pol(Pol *pol, double start, double epsilon){
	// Returns a root of the polynomial pol closed to start, found with the Haley's method with a
	// precision epsilon, in the interval ]0,1[.

	Pol dpol; // first derivative of pol
	init_pol(&dpol,1);
	der_pol(&dpol,pol);

	Pol ddpol; // second derivative of pol
	init_pol(&ddpol,1);
	der_pol(&ddpol,&dpol);

	double dp; // pol'(x)
	double pdp; // pol(x) / pol'(x)
	while (pow(eval_pol(pol,start),2) > pow(epsilon,2) || start <= 0 || start >= 1){
		dp = eval_pol(&dpol,start);
		pdp = eval_pol(pol,start)/dp;
		start -= pdp / (1 - pdp * eval_pol(&ddpol,start) / (2 * dp));
	}

	free(dpol.coef);
	free(ddpol.coef);

	return start;
}
