#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "init.h"
#include "param.h"
#include "ellipsoids.h"
#include "maths.h"
#include "integration.h"
#include "measure.h"

// int main(){

// 	// matrix 3x3 product

// 	// double mat1[9] = {0,1,0,1,0,0,0,0,1};
// 	// double mat2[9] = {0,1,0,1,0,0,0,0,1};

// 	// Pol pol1[9], pol2[9];
// 	// init_pol(&pol1[0],9);
// 	// init_pol(&pol2[0],9);

// 	// scal2pol(&pol1[0],&mat1[0],9);
// 	// scal2pol(&pol2[0],&mat2[0],9);

// 	// Pol pol[9];
// 	// //init_pol(&pol[0],4);

// 	// mat_pol_prod(&pol[0],&pol1[0],&pol2[0],3,3,3);
// 	// for (int i = 0; i < 3; i++){
// 	// 	for (int j = 0; j < 3; j++){
// 	// 		printf("%f ",pol[i * 3 + j].coef[0]);
// 	// 	}
// 	// 	printf("\n");
// 	// }

// 	// for (int c = 0; c < 9; c++){
// 	// 	free(pol[c].coef);
// 	// 	free(pol1[c].coef);
// 	// 	free(pol2[c].coef);
// 	// }

// 	// matrix 3x3 determinant

// 	// double mat1[9] = {2,2,1,0,3,1,-1,3,4};
// 	// Pol pol,mat[9];

// 	// scal2pol(&mat[0],&mat1[0],9);

// 	// det3_mat_pol(&pol,&mat[0]);
// 	// printf("%f\n",pol.coef[0]);

// 	// free(pol.coef);
// 	// for (int i = 0; i < 9; i++){
// 	// 	free(mat[i].coef);
// 	// }

// 	// characteristic polynomial matrix 2x2

// 	// Pol pol[4];
// 	// double a[2] = {-1,1};
// 	// create_pol(&pol[0],1,&a[0]);
// 	// double b[1] = {0};
// 	// create_pol(&pol[1],0,&b[0]);
// 	// double c[1] = {0};
// 	// create_pol(&pol[2],0,&c[0]);
// 	// double d[2] = {-1,1};
// 	// create_pol(&pol[3],1,&d[0]);

// 	// Pol det;
// 	// det2_mat_pol(&det,&pol[0]);

// 	// for (int c = 0; c <= det.deg; c++){
// 	// 	printf("%f ",det.coef[c]);
// 	// }
// 	// printf("\n");

// 	// free_mat_pol(&pol[0],4);
// 	// free(det.coef);

// 	// adjoint of 3*3 matrix

// 	// double mat1[9] = {2,2,1,0,3,1,-1,3,4};
// 	// Pol adj[9],mat[9];

// 	// scal2pol(&mat[0],&mat1[0],9);

// 	// adj_mat3_pol(&adj[0],&mat[0]);
// 	// for (int i = 0; i < 3; i++){
// 	// 	for (int j = 0; j < 3; j++){
// 	// 		printf("%f ",adj[i*3+j].coef[0]);
// 	// 	}
// 	// 	printf("\n");
// 	// }

// 	// free_mat_pol(&adj[0],9);
// 	// free_mat_pol(&mat[0],9);
// }

// int main(){

// 	Ellipsoid el1,el2;
// 	SizeEllipsoid sizes[2];

// 	double ax1[3] = {1,1,1};
// 	double ax2[3] = {1,1,1};
// 	equal(sizes[0].ax,ax1,3);
// 	equal(sizes[1].ax,ax2,3);
// 	double ax1inv2[3];
// 	double ax2inv2[3];
// 	for (int coord = 0; coord < 3; coord++){
// 		ax1inv2[coord] = 1/pow(ax1[coord],2);
// 		ax2inv2[coord] = 1/pow(ax2[coord],2);
// 	}
// 	sizes[0].g = 1/ax1inv2[max(ax1inv2,3)];
// 	sizes[1].g = 1/ax2inv2[max(ax2inv2,3)];

// 	double r1[3] = {0,0,0};
// 	double r2[3] = {10,0,0};
// 	double v[3] = {0,0,0};
// 	double q1[4] = {0,0,0,1};
// 	double q2[4] = {0,0,0,1};

// 	equal(el1.r,r1,3);
// 	equal(el2.r,r2,3);
// 	equal(el1.v,v,3);
// 	equal(el2.v,v,3);
// 	equal(el1.q,q1,4);
// 	equal(el2.q,q2,4);
// 	el1.size = 0;
// 	el2.size = 1;

// 	init(el1.Q,9);
// 	init(el2.Q,9);
// 	rotation_matrix(&el1);
// 	rotation_matrix(&el2);

// 	printf("r1: %f %f %f\n",r1[0],r1[1],r1[2]);
// 	printf("r2: %f %f %f\n",r2[0],r2[1],r2[2]);
// 	printf("q1: %f %f %f %f\n",q1[0],q1[1],q1[2],q1[3]);
// 	printf("q2: %f %f %f %f\n",q2[0],q2[1],q2[2],q2[3]);

// 	double mu = scaling(&el1,&el2,&sizes[0]);
// 	printf("mu: %f\n",mu);
// }

int main(int argc, char *argv[]){

	double mat[9];
	double mat1[9] = {1,2,-1,4,2,-2,1,1,2};
	double mat2[9] = {-1,-2,-3,2,-1,1,1,3,3};
	double mat3[3] = {1,-2,4};

	matrix_prod_333(&mat[0],&mat1[0],&mat2[0]);
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			printf("%f ",mat[i*3+j]);
		}
		printf("\n");
	}
	printf("\n");

	matrix_prod_313(&mat[0],&mat3[0],&mat3[0]);
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			printf("%f ",mat[i*3+j]);
		}
		printf("\n");
	}
	printf("\n");

	matrix_prod_331(&mat[0],&mat1[0],&mat3[0]);
	for (int i = 0; i < 3; i++){
		printf("%f ",mat[i]);
	}
	printf("\n");

	// Ellipsoid el1,el2;
	// SizeEllipsoid sizes[2];
	// Overlap overlap;
	// Par par;
	// Forces forces[2];
	// for (int el = 0; el < 2; el++){
	// 	for (int coord = 0; coord < 3; coord++){
	// 		forces[el].el[coord] = 0;
	// 		forces[el].dis[coord] = 0;
	// 		forces[el].M[coord] = 0;
	// 	}
	// }
	// par.epsilon_dist2 = 1e-8;
	// par.epsilon_ang = 1e-3;
	// par.ke = 1;
	// par.kd = 1;

	// // double ax1[3] = {0.5,0.5,0.5};
	// // double ax2[3] = {0.5,0.5,0.5};
	// double ax1[3] = {1,0.5,0.5};
	// double ax2[3] = {1,0.5,0.5};
	// equal(sizes[0].ax,ax1,3);
	// equal(sizes[1].ax,ax2,3);
	// double ax1inv2[3];
	// double ax2inv2[3];
	// for (int coord = 0; coord < 3; coord++){
	// 	ax1inv2[coord] = 1/pow(ax1[coord],2);
	// 	ax2inv2[coord] = 1/pow(ax2[coord],2);
	// }
	// sizes[0].g = 1/ax1inv2[max(ax1inv2,3)];
	// sizes[1].g = 1/ax2inv2[max(ax2inv2,3)];

	// // double r1[3] = {-1.44402,1.69752,0.845323};
	// // double r2[3] = {-2.0953,1.02594,1.17792};
	// double r1[3] = {0.1,0,0};
	// double r2[3] = {0.8,0.4,0};
	// double v[3] = {0,0,0};
	// double q1[4] = {0,0,0,1};
	// double q2[4] = {0,0,-0.3090169943749474,0.9510565162951535};
	// double qdot[4] = {0,0,0,0};
	// int s = 0;

	// equal(el1.r,r1,3);
	// equal(el2.r,r2,3);
	// equal(el1.v,v,3);
	// equal(el2.v,v,3);
	// equal(el1.q,q1,4);
	// equal(el2.q,q2,4);
	// // equal(el1.qdot,qdot,4);
	// // equal(el2.qdot,qdot,4);
	// el1.size = 0;
	// el2.size = 1;

	// init(el1.Q,9);
	// init(el2.Q,9);
	// rotation_matrix(&el1);
	// rotation_matrix(&el2);

	// printf("r1: %f %f %f\n",r1[0],r1[1],r1[2]);
	// printf("r2: %f %f %f\n",r2[0],r2[1],r2[2]);
	// // printf("q1: %f %f %f %f\n",q1[0],q1[1],q1[2],q1[3]);
	// // printf("q2: %f %f %f %f\n",q2[0],q2[1],q2[2],q2[3]);

	// equal(overlap.c1,r1,3);
	// equal(overlap.c2,r2,3);
	// overlap.mu = 1;
	// par.epsilon_mu = 1;
	// par.epsilon = 1e-7;
	// contact(&par,&overlap,&sizes[0],&el1,&el2);
	// double n[3];
	// n_surface_vec(&n[0],&el1,&sizes[0],&overlap.r_c[0]);
	// printf("n: %f %f %f\n",n[0],n[1],n[2]);
	// // for (int coord = 0; coord < 3; coord++){
	// // 	overlap.r_c[coord] = el1.r[coord] + 0.5*(el2.r[coord] - el1.r[coord]);
	// // }
	// // all_forces(&par,&forces[0],&forces[1],&el1,&el2,&sizes[0],&overlap);

	// // double diff[3];
	// // matrix_diff(&diff[0],&el2.r[0],&el1.r[0],3);
	// // double f = 1 - matrix_norm(&diff[0],3);

	// printf("overlaping: %d mu: %f\nr_c: %f %f %f\nr_c1: %f %f %f\nr_c2: %f %f %f\n",overlap.overlapping,overlap.mu,overlap.r_c[0],overlap.r_c[1],overlap.r_c[2],overlap.r_c1[0],overlap.r_c1[1],overlap.r_c1[2],overlap.r_c2[0],overlap.r_c2[1],overlap.r_c2[2]);
	// // printf("fel12: %f %f %f\n",forces[0].el[0],forces[0].el[1],forces[0].el[2]);
	// // printf("fth: %f %f %f\n",f*diff[0],f*diff[1],f*diff[2]);
	// //printf("%f\n",belonging_function(&el,r));
	// // scaling(&par,&el1,&el2,&sizes[0],&overlap);
	// // printf("NEW METHOD:\noverlaping: %d mu: %f\nr_c: %f %f %f\nr_c1: %f %f %f\nr_c2: %f %f %f\n",overlap.overlapping,overlap.mu,overlap.r_c[0],overlap.r_c[1],overlap.r_c[2],overlap.r_c1[0],overlap.r_c1[1],overlap.r_c1[2],overlap.r_c2[0],overlap.r_c2[1],overlap.r_c2[2]);

}
