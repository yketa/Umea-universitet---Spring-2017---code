#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "measure.h"

FILE* init_measure(int argc, char *argv[], Par *par, SizeEllipsoid *sizes, Ellipsoid *ellipsoids, Forces *forces,
	Overlap *overlaps, double *p, double *L){
	// Initialises the meausrement file.

	FILE* file = NULL; // results.csv file

	if (argc == 1){
		// no name given by user
    	file = fopen("results.csv", "w");
    }
    else{
    	// name given by user
    	file = fopen(argv[1],"w");
    }

    if (file == NULL){
        // error opening file
        return file;
    }

	fprintf(file,"time,"); // time

	for (int el = 0; el < nEllipsoids; el++){
		fprintf(file,"x%d,y%d,z%d,",el,el,el); // position
		fprintf(file,"vx%d,vy%d,vz%d,",el,el,el); // velocity
		fprintf(file,"wx%d,wy%d,wz%d,",el,el,el); // angular velocity
		fprintf(file,"felx%d,fely%d,felz%d,",el,el,el); // elastic force
		fprintf(file,"fdisx%d,fdisy%d,fdisz%d,",el,el,el); // dissipative force
		fprintf(file,"orx%d,ory%d,orz%d,",el,el,el); // orientation of the longest axis
	}

	int couple = 0;
	for (int el1 = 0; el1 < nEllipsoids - 1; el1++){
		for (int el2 = el1 + 1; el2 < nEllipsoids; el2++){
			fprintf(file,"v%dc%dx,v%dc%dy,v%dc%dz,",el1,el2,el1,el2,el1,el2); // velocity of particle el1 at its point of contact with particle el2
			fprintf(file,"v%dc%dx,v%dc%dy,v%dc%dz,",el2,el1,el2,el1,el2,el1); // velocity of particle el2 at its point of contact with particle el1
			fprintf(file,"rc%d_in%dx,rc%d_in%dy,rc%d_in%dz,",el1,el1,el1,el1,el1,el1); // relative squared coordinates of the contact point on particle el1
			fprintf(file,"rc%d_in%dx,rc%d_in%dy,rc%d_in%dz,",el2,el2,el2,el2,el2,el2); // relative squared coordinates of the contact point on particle el2
			fprintf(file,"err%d%d,",el1,el2); // error on the determination of the root of polynomial H_{el1el2}
			fprintf(file,"overlap%d%d,",el1,el2); // overlap between el1 and el2
			couple++;
		}
	}

	fprintf(file,"dpx,dpy,dpz,dLx,dLy,dLz,"); // variations of momentum and angular momentum

	fprintf(file,"\n");

	par->iter = -1;
	measure(file, par, sizes, ellipsoids, forces, overlaps, p, L); // values at t = 0
		// note that the first elements dp and dL will correspond to the initial values of
		// the total translational momentum and translational momentum

	return file;

}

void measure(FILE* file, Par *par, SizeEllipsoid *sizes, Ellipsoid *ellipsoids, Forces *forces, Overlap *overlaps,
	double *p, double *L){
	// Prints in .csv format the time and the positions, velocities and angular velocities of the ellipsoids.

	double time = (par->iter+1)*dt; // time
	fprintf(file,"%f,",time); // printing time

	double ax[3]; // direction of the longest axis

	for (int el = 0; el < nEllipsoids; el++){

		init(&ax[0],3);
		ax[(sizes + (ellipsoids + el)->size)->axMax] = 1; // longest axis
		action_q(&ax[0],(ellipsoids + el)->q,&ax[0]); // rotated longest axis

		fprintf(file,"%f,%f,%f,",(ellipsoids + el)->r[0],(ellipsoids + el)->r[1],(ellipsoids + el)->r[2]); // position
		fprintf(file,"%f,%f,%f,",(ellipsoids + el)->v[0],(ellipsoids + el)->v[1],(ellipsoids + el)->v[2]); // velocity
		fprintf(file,"%f,%f,%f,",(ellipsoids + el)->w[0],(ellipsoids + el)->w[1],(ellipsoids + el)->w[2]); // angular velocity
		fprintf(file,"%f,%f,%f,",(forces + el)->el[0],(forces + el)->el[1],(forces + el)->el[2]); // elastic force
		fprintf(file,"%f,%f,%f,",(forces + el)->dis[0],(forces + el)->dis[1],(forces + el)->dis[2]); // dissipative force
		fprintf(file,"%f,%f,%f,",ax[0],ax[1],ax[2]); // orientation of the longest axis
	}

	int couple = 0;
	for (int el1 = 0; el1 < nEllipsoids - 1; el1++){
		for (int el2 = el1 + 1; el2 < nEllipsoids; el2++){
			if ((overlaps + couple)->overlapping){
				// particles are overlapping
				fprintf(file,"%f,%f,%f,",(overlaps + couple)->v1c2[0],(overlaps + couple)->v1c2[1],(overlaps + couple)->v1c2[2]); // velocity of particle el1 at its point of contact with particle el2
				fprintf(file,"%f,%f,%f,",(overlaps + couple)->v2c1[0],(overlaps + couple)->v2c1[1],(overlaps + couple)->v2c1[2]); // velocity of particle el2 at its point of contact with particle el1
				fprintf(file,"%f,%f,%f,",(overlaps + couple)->r_c1_in1[0],(overlaps + couple)->r_c1_in1[1],(overlaps + couple)->r_c1_in1[2]); // relative squared coordinates of the contact point on particle el1
				fprintf(file,"%f,%f,%f,",(overlaps + couple)->r_c2_in2[0],(overlaps + couple)->r_c2_in2[1],(overlaps + couple)->r_c2_in2[2]); // relative squared coordinates of the contact point on particle el2
				fprintf(file,"%f,",(overlaps + couple)->err); // error on the determination of the root of polynomial H_{el1el2}
				fprintf(file,"%f,",1 - (overlaps + couple)->mu); // overlap between el1 and el2
			}
			else{
				// particles are not overlapping
				fprintf(file,"%d,%d,%d,",0,0,0); // velocity of particle el1 at its point of contact with particle el2
				fprintf(file,"%d,%d,%d,",0,0,0); // velocity of particle el2 at its point of contact with particle el1
				fprintf(file,"%d,%d,%d,",0,0,0); // relative squared coordinates of the contact point on particle el1
				fprintf(file,"%d,%d,%d,",0,0,0); // relative squared coordinates of the contact point on particle el2
				fprintf(file,"%d,",0); // error on the determination of the root of polynomial H_{el1el2}
				fprintf(file,"%d,",0); // overlap between el1 and el2
			}
			couple++;
		}
	}

	momentum(ellipsoids,sizes,p,L); // momentum and angular momentum

	for (int coord = 0; coord < 3; coord++){
		// variation of momentum
		fprintf(file,"%f,",p[coord + 3] - p[coord]);
	}
	for (int coord = 0; coord < 3; coord++){
		// variation of angular momentum
		fprintf(file,"%f,",L[coord + 3] - L[coord]);
	}

	fprintf(file,"\n");

}

void momentum(Ellipsoid *ellipsoids, SizeEllipsoid *sizes, double *p, double *L){
	// Associates to p and L the total translational and angular momentum on each axis
	// on two consecutives time steps.

	for (int coord = 0; coord < 3; coord++){
		// passing to new iteration
		p[coord] = p[coord + 3];
		p[coord + 3] = 0;
		L[coord] = L[coord + 3];
		L[coord + 3] = 0;
	}

	double T[3], TR[3], Rb[3], R[3], A[3]; // translational momentum, translational angular momentum, rotation angular momentum in the body and world frame and total angular momentum

	for (int el = 0; el < nEllipsoids; el++){

		scal_prod(&T[0],(sizes + (ellipsoids + el)->size)->m,(ellipsoids + el)->v,3); // translational momentum of ellipsoid el

		for (int coord = 0; coord < 3; coord++){
			Rb[coord] = (sizes + (ellipsoids + el)->size)->I[coord] * (ellipsoids + el)->w[coord]; // rotational angular momentum of ellipsoid el in the body frame
		}

		for (int i = 0; i < 3; i++){
			R[i] = 0;
			for (int j = 0; j < 3; j++){
				R[i] += (ellipsoids + el)->Q[i*3 + j]*Rb[j]; // rotational angular momentum of ellipsoid el in the world frame
			}
		}

		cross(&TR[0],(ellipsoids + el)->r,&T[0]); // translational angular momentum
		matrix_sum(&A[0],&TR[0],&R[0],3); // total angular momentum

		for (int coord = 0; coord < 3; coord++){
			p[coord + 3] += T[coord]; // total translational momentum
			L[coord + 3] += A[coord]; // total angular momentum
		}

	}

}
