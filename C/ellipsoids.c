#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "ellipsoids.h"

void rotation_matrix(Ellipsoid *ellipsoid){
	// Updates the rotation matrix associated to ellipsoid->q.

	rotation(ellipsoid->Q,ellipsoid->q); // rotation matrix associated to ellipsoid->q
}

double belonging_function(Ellipsoid *ellipsoid, SizeEllipsoid *sizes, double *vec){
	// Returns the belonging fonction of ellipsoid, with semi-axes ax, evaluated in vec.

	double *ax = &(sizes + ellipsoid->size)->ax[0]; // semi-axes of the ellipsoid

	double rot_vec[3];

	double result = 0; // evaluation of the belonging function in vec
	for (int i = 0; i < 3; i++){
		rot_vec[i] = 0;
		for (int j = 0; j < 3; j++){
			rot_vec[i] += (vec[j] - ellipsoid->r[j]) * ellipsoid->Q[j*3 + i]; // = ellipsoid->Q^T (vec - ellipsoid->r)
		}
		result += pow(rot_vec[i],2)/pow(ax[i],2); // = ellipsoid->Q diag(R_i^{-2}) ellipsoid->Q^T (vec - ellipsoid->r)
	}

	result -= 1;
	return result;
}

void qdot(double *qp, Ellipsoid *ellipsoid){
	// Associates qp to the first derivative of the quaternion of ellipsoid.

	double wquat[4]; // quaternion associated to w
	equal(&wquat[0],ellipsoid->w,3);
	wquat[3] = 0;

	qprod(qp,ellipsoid->q,&wquat[0]);
	scal_prod(qp,0.5,qp,4);
}

void wdot(Ellipsoid *ellipsoid, SizeEllipsoid *sizes, Forces *forces){
	// Associates to forces->wdot the first derivative of the rotation vector of the ellipsoid.

	double *I = (sizes + ellipsoid->size)->I; // moments of inertia of ellipsoid

	double Iw[3]; // product I w
	for (int coord = 0; coord < 3; coord++){
		Iw[coord] = I[coord]*ellipsoid->w[coord];
	}

	double wIw[3]; // cross product w x Iw
	cross(&wIw[0],ellipsoid->w,&Iw[0]);

	for (int coord = 0; coord < 3; coord++){
		forces->wdot[coord] = (forces->M[coord] - wIw[coord])/I[coord]; // \dot{w}
	}

}

void n_surface_vec(double *n, Ellipsoid *ellipsoid, SizeEllipsoid *sizes, double *vec){
	// Associates n to the non-unit surface vector of ellipsoid, whose semi-axes are ax,
	// at vec.

	double *ax = &(sizes + ellipsoid->size)->ax[0]; // semi-axes of the ellipsoid

	double n_[3];
	for (int i = 0; i < 3; i++){
		n_[i] = 0;
		for (int j = 0; j < 3; j++){
			n_[i] += (vec[j] - ellipsoid->r[j]) * ellipsoid->Q[j*3 + i]; // = ellipsoid->Q^T (vec - ellipsoid1->r)
		}
		n_[i] /= pow(ax[i],2); // = diag(1/ax^2) ellipsoid->Q^T (vec - ellipsoid->r)
	}

	matrix_prod_331(&n[0],ellipsoid->Q,&n_[0]); // = \underbrace{ellipsoid->Q diag(1/ax^2) ellipsoid->Q^T}_{\bar{B}} (vec - ellipsoid->r)
}

void inv_red_bel(double *Binv, Ellipsoid *ellipsoid, SizeEllipsoid *sizes){
	// Associates Binv to the inverse of the reduced belonging matrix of ellipsoid.

	double mat9[9]; // useful matrix

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			mat9[j * 3 + i] = ellipsoid->Q[i * 3 + j]*pow((sizes + ellipsoid->size)->ax[j],2); // diag(R^2) tQ
		}
	}
	matrix_prod_333(&Binv[0],ellipsoid->Q,&mat9[0]); // Q diag(R^2) tQ = Binv
}

void adj_y_ab(double *adj, double *Bainv, double *Bbinv){
	// Associates to adj the adjoint of the matrix Y_{AB} defined by the inverses of
	// the reduced belonging matrix of the ellipsoids, Bainv and Bbinv.
	// Coefficients calculated with Mathematica.

   adj[0] = - Bainv[5]*Bainv[7] + Bainv[4]*Bainv[8];
   adj[1] = Bainv[8]*Bbinv[4] - Bainv[7]*Bbinv[5] - Bainv[5]*Bbinv[7] + Bainv[4]*Bbinv[8];
   adj[2] = adj[0] - adj[1] - Bbinv[5]*Bbinv[7] + Bbinv[4]*Bbinv[8];
   adj[1] -= 2*adj[0];

   adj[3] = Bainv[2]*Bainv[7] - Bainv[1]*Bainv[8];
   adj[4] = - Bainv[8]*Bbinv[1] + Bainv[7]*Bbinv[2] + Bainv[2]*Bbinv[7] - Bainv[1]*Bbinv[8];
   adj[5] = adj[3] - adj[4] + Bbinv[2]*Bbinv[7] - Bbinv[1]*Bbinv[8];
   adj[4] -= 2*adj[3];

   adj[6] = - Bainv[2]*Bainv[4] + Bainv[1]*Bainv[5];
   adj[7] = Bainv[5]*Bbinv[1] - Bainv[4]*Bbinv[2] - Bainv[2]*Bbinv[4] + Bainv[1]*Bbinv[5];
   adj[8] = adj[6] - adj[7] - Bbinv[2]*Bbinv[4] + Bbinv[1]*Bbinv[5];
   adj[7] -= 2*adj[6];

   adj[9] = Bainv[5]*Bainv[6] - Bainv[3]*Bainv[8];
   adj[10] = - Bainv[8]*Bbinv[3] + Bainv[6]*Bbinv[5] + Bainv[5]*Bbinv[6] - Bainv[3]*Bbinv[8];
   adj[11] = adj[9] - adj[10] + Bbinv[5]*Bbinv[6] - Bbinv[3]*Bbinv[8];
   adj[10] -= 2*adj[9];

   adj[12] = - Bainv[2]*Bainv[6] + Bainv[0]*Bainv[8];
   adj[13] = Bainv[8]*Bbinv[0] - Bainv[6]*Bbinv[2] - Bainv[2]*Bbinv[6] + Bainv[0]*Bbinv[8];
   adj[14] = adj[12] - adj[13] - Bbinv[2]*Bbinv[6] + Bbinv[0]*Bbinv[8];
   adj[13] -= 2*adj[12];

   adj[15] = Bainv[2]*Bainv[3] - Bainv[0]*Bainv[5];
   adj[16] = - Bainv[5]*Bbinv[0] + Bainv[3]*Bbinv[2] + Bainv[2]*Bbinv[3] - Bainv[0]*Bbinv[5];
   adj[17] = adj[15] - adj[16] + Bbinv[2]*Bbinv[3] - Bbinv[0]*Bbinv[5];
   adj[16] -= 2*adj[15];

   adj[18] = - Bainv[4]*Bainv[6] + Bainv[3]*Bainv[7];
   adj[19] = Bainv[7]*Bbinv[3] - Bainv[6]*Bbinv[4] - Bainv[4]*Bbinv[6] + Bainv[3]*Bbinv[7];
   adj[20] = adj[18] - adj[19] - Bbinv[4]*Bbinv[6] + Bbinv[3]*Bbinv[7];
   adj[19] -= 2*adj[18];

   adj[21] = Bainv[1]*Bainv[6] - Bainv[0]*Bainv[7];
   adj[22] = - Bainv[7]*Bbinv[0] + Bainv[6]*Bbinv[1] + Bainv[1]*Bbinv[6] - Bainv[0]*Bbinv[7];
   adj[23] = adj[21] - adj[22] + Bbinv[1]*Bbinv[6] - Bbinv[0]*Bbinv[7];
   adj[22] -= 2*adj[21];

   adj[24] = -Bainv[1]*Bainv[3] + Bainv[0]*Bainv[4];
   adj[25] = Bainv[4]*Bbinv[0] - Bainv[3]*Bbinv[1] - Bainv[1]*Bbinv[3] + Bainv[0]*Bbinv[4];
   adj[26] = adj[24] - adj[25] - Bbinv[1]*Bbinv[3] + Bbinv[0]*Bbinv[4];
   adj[25] -= 2*adj[24];

}

void p_ab(double *pab, double *adj, double *rab){
	// Associates to pab the coefficients of the polynomial P_{AB} corresponding to
	// the ajoint adj of Y_{AB} and ellipsoids whose centers are separated by vector
	// rab = \vec{r_B} - \vec{r_A}.

	// DETERMINATION OF THE POLYNOMIAL P_{AB}

	init(pab,5); // initialisation of the coefficients of the polynomial

	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			for (int co = 0; co <= 2; co++){
				pab[co + 1] += adj[(i * 3 + j)*3 + co]*rab[i]*rab[j]; // lambda (rB - rA)^T adj (rB - rA)
				pab[co + 2] -= adj[(i * 3 + j)*3 + co]*rab[i]*rab[j]; // -lambda^2 (rB - rA)^T adj (rB - rA)
			}
		}
	}

}

void q_ab(double *det, double *Bainv, double *Bbinv){
	// Associated to det the coefficients of the determinant of the matrix Y_{AB}
	// defined the inverses of the reduced belonging matrix of the ellipsoids, Bainv
	// and Bbinv.
	// Coefficients calculated with Mathematica.

	// USEFUL ITEMS

	double aaa = - Bainv[2]*Bainv[4]*Bainv[6] + Bainv[1]*Bainv[5]*Bainv[6] + Bainv[2]*Bainv[3]*Bainv[7] - Bainv[0]*Bainv[5]*Bainv[7] - Bainv[1]*Bainv[3]*Bainv[8] + Bainv[0]*Bainv[4]*Bainv[8];

	double aab = - Bainv[5]*Bainv[7]*Bbinv[0] + Bainv[4]*Bainv[8]*Bbinv[0] + Bainv[5]*Bainv[6]*Bbinv[1] - Bainv[3]*Bainv[8]*Bbinv[1] - Bainv[4]*Bainv[6]*Bbinv[2] + Bainv[3]*Bainv[7]*Bbinv[2];
	aab += Bainv[2]*Bainv[7]*Bbinv[3] - Bainv[1]*Bainv[8]*Bbinv[3] - Bainv[2]*Bainv[6]*Bbinv[4] + Bainv[0]*Bainv[8]*Bbinv[4] + Bainv[1]*Bainv[6]*Bbinv[5] - Bainv[0]*Bainv[7]*Bbinv[5];
	aab += - Bainv[2]*Bainv[4]*Bbinv[6] + Bainv[1]*Bainv[5]*Bbinv[6] + Bainv[2]*Bainv[3]*Bbinv[7] - Bainv[0]*Bainv[5]*Bbinv[7] - Bainv[1]*Bainv[3]*Bbinv[8] + Bainv[0]*Bainv[4]*Bbinv[8];

	double abb = - Bainv[8]*Bbinv[1]*Bbinv[3] + Bainv[7]*Bbinv[2]*Bbinv[3]+ Bainv[8]*Bbinv[0]*Bbinv[4] - Bainv[6]*Bbinv[2]*Bbinv[4] - Bainv[7]*Bbinv[0]*Bbinv[5] + Bainv[6]*Bbinv[1]*Bbinv[5];
	abb += Bainv[5]*Bbinv[1]*Bbinv[6] - Bainv[4]*Bbinv[2]*Bbinv[6] - Bainv[2]*Bbinv[4]*Bbinv[6] + Bainv[1]*Bbinv[5]*Bbinv[6] - Bainv[5]*Bbinv[0]*Bbinv[7] + Bainv[3]*Bbinv[2]*Bbinv[7];
	abb += Bainv[2]*Bbinv[3]*Bbinv[7] - Bainv[0]*Bbinv[5]*Bbinv[7] + Bainv[4]*Bbinv[0]*Bbinv[8] - Bainv[3]*Bbinv[1]*Bbinv[8] - Bainv[1]*Bbinv[3]*Bbinv[8] + Bainv[0]*Bbinv[4]*Bbinv[8];

	double bbb = - Bbinv[2]*Bbinv[4]*Bbinv[6] + Bbinv[1]*Bbinv[5]*Bbinv[6] - Bbinv[0]*Bbinv[5]*Bbinv[7] + Bbinv[2]*Bbinv[3]*Bbinv[7] - Bbinv[1]*Bbinv[3]*Bbinv[8] + Bbinv[0]*Bbinv[4]*Bbinv[8];

	// DETERMINATION OF POLYNOMIAL Q_{AB}

	det[0] = aaa;
	det[1] = -3*aaa + aab;
	det[2] = 3*aaa - 2*aab + abb;
	det[3] = - aaa + aab - abb + bbb;

}

double root_Haley_pol_hAB(double *hab, double *dhab, double start, double epsilon){
	// Returns the root of the polynomial H_{AB} closed to start, found with the Haley's method with a
	// precision epsilon, in the interval ]0,1[.

	double ddhab[5]; // second derivative of H_{AB}
	for (int co = 0; co<= 4; co++){
		ddhab[co] = (co + 1)*(co + 2)*hab[co + 2]; // coefficients of the polynomial
	}

	double p = eval_pol(hab,6,start); // H_{AB}(x)
	double dp; // H_{AB}'(x)
	while (p < - epsilon || p > epsilon || start <= 0 || start >= 1){
		dp = eval_pol(dhab,5,start);
		start -= (p/dp) / (1 - (p/dp) * eval_pol(&ddhab[0],4,start) / (2 * dp));
		p = eval_pol(hab,6,start);
	}

	return start;
}

double eval_pol(double *coef, int deg, double x){
	// Evaluates the polynomial of coefficients coef, degree deg
	// in x using the Horner's rule.

	double result = coef[deg]; // a_n
	for (int c = deg - 1; c >= 0; c--){
		result = coef[c] + result * x; // a_{n-1} + b_n x
	}

	return result;
}

void contact(Par *par, Overlap *overlap, SizeEllipsoid *sizes,
	Ellipsoid *ellipsoid1, Ellipsoid *ellipsoid2){
	// Updates the contact points between ellipsoid1 and ellipsoid2.
	// Algorithm described in subsection 2.5.1 — based on Perram and Wertheim, J. of Comp. Phys., 1995

	// DETERMINATION OF THE INVERSE OF THE REDUCED BELONGING MATRIX

	double Bainv[9]; // inverse of the reduced belong matrix of ellipsoid1
	inv_red_bel(&Bainv[0],ellipsoid1,sizes);

	double Bbinv[9]; // inverse of the reduced belong matrix of ellipsoid2
	inv_red_bel(&Bbinv[0],ellipsoid2,sizes);

	// DETERMINATION OF THE ADJOINT OF POLYNOMIAL Y_{AB}

	double adj[27]; // coefficients of polynomial Y_{AB}
	adj_y_ab(&adj[0],&Bainv[0],&Bbinv[0]);

	// DETERMINATION OF THE DETERMINANT OF Y_{AB}

	double qab[4]; // coefficients of the determinant of matrix Y_{AB} (polynomial Q_{AB})
	q_ab(&qab[0],&Bainv[0],&Bbinv[0]);

	// DETERMINATION OF R_{AB}

	double rab[3]; // \vec{r_B} - \vec{r_A}
	matrix_diff(&rab[0],ellipsoid2->r,ellipsoid1->r,3);

	// DETERMINATION OF POLYNOMIAL P_{AB}

	double pab[5]; // coefficients of the polynomial P_{AB}
	p_ab(&pab[0],&adj[0],&rab[0]);

	// DETERMINATION OF POLYNOMIALS H_{AB} AND H'_{AB}

	double hab[7]; // coefficients of the polynomial H_{AB}
	init(&hab[0],7); // initialisation of the coefficients of the polynomial

	for (int copp = 1; copp <= 4; copp++){
		for (int coq = 0; coq <= 3; coq++){
			hab[copp + coq - 1] += copp*pab[copp]*qab[coq]; // P'_{AB} Q_{AB}
		}
	}

	for (int cop = 0; cop <= 4; cop++){
		for (int coqq = 1; coqq <= 3; coqq++){
			hab[cop + coqq - 1] -= pab[cop]*coqq*qab[coqq]; // - P_{AB} Q'_{AB}
		}
	}

	double dhab[6]; // coefficients of the polynomial H'_{AB}
	
	for (int co = 0; co <= 5; co++){
		dhab[co] = (co + 1)*hab[co + 1]; // coefficients of the derivative of hab
	}

	// DETERMINATION OF THE ROOT OF H_{AB} IN [0,1]

	double root = (sizes + ellipsoid1->size)->R/((sizes + ellipsoid1->size)->R + (sizes + ellipsoid2->size)->R); // initial guess for the root of H_{AB}
	root = root_Haley_pol_hAB(&hab[0],&dhab[0],root,par->epsilon); // root in the interval [0,1]

	// ERROR CHECKING

	overlap->err = 2 * par->epsilon / eval_pol(&dhab[0],5,root); // length of the interval where the root could have been found

	// DETERMINATION OF THE RESCALING FACTOR

	overlap->mu = sqrt(eval_pol(&pab[0],4,root)/eval_pol(&qab[0],3,root)); // mu = sqrt(potential)

	if (overlap->mu < 1){
		// ellipsoids are overlapping
		overlap->overlapping = 1;

		// CONTACT POINTS

		double det = eval_pol(&qab[0],3,root); // determinant of Y_{AB}(root)

		double yabinvrab[3]; // Y_{AB}^-1 (\vec{r_B} - \vec{r_A})
		init(&yabinvrab[0],3); // initialisation
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				yabinvrab[i] += eval_pol(&adj[(i * 3 + j)*3],2,root) * rab[j] / det; // adj(Y_{AB}) (\vec{r_B} - \vec{r_A}) / det(Y_{AB}) = Y_{AB}^-1 (\vec{r_B} - \vec{r_A})
			}
		}

		init(overlap->r_c1_in1,3); // initialisation of the relative squared coordinates of the contact point on particle 1
		init(overlap->r_c2_in2,3); // initialisation of the relative squared coordinates of the contact point on particle 2

		for (int i = 0; i < 3; i++){
			// contact point
			overlap->r_c[i] = ellipsoid2->r[i]; // r_B
			for (int j = 0; j < 3; j++){
				overlap->r_c[i] -= root * Bbinv[i * 3 + j] * yabinvrab[j]; // - lambda B_B^-1 Y_{AB}^-1 (\vec{r_B} - \vec{r_A})
			}
			// contact point on particles
			overlap->r_c1[i] = ellipsoid1->r[i] + (overlap->r_c[i] - ellipsoid1->r[i])/overlap->mu; // contact point on particle 1
			overlap->r_c2[i] = ellipsoid2->r[i] + (overlap->r_c[i] - ellipsoid2->r[i])/overlap->mu; // contact point on particle 2
			// relative squared coordinates of the point on particles
			for (int l = 0; l < 3; l++){
				overlap->r_c1_in1[l] += ellipsoid1->Q[i * 3 + l] * (overlap->r_c1[i] - ellipsoid1->r[i]); // tQa (r_ca - ra)
				overlap->r_c2_in2[l] += ellipsoid2->Q[i * 3 + l] * (overlap->r_c2[i] - ellipsoid2->r[i]); // tQb (r_cb - rb)
			}
		}

	}
	else {
		// ellipsoids are not overlapping
		overlap->overlapping = 0;
	}

}
