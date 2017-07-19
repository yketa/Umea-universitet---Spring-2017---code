#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "init.h"
#include "param.h"
#include "ellipsoids.h"
#include "maths.h"
#include "integration.h"
#include "measure.h"

int iter; // iteration number

int main(int argc, char *argv[]){

	// PARTICLES

	Ellipsoid ellipsoids[nEllipsoids]; // ellipsoids
	Forces forces[nEllipsoids]; // forces exerted on the ellipsoids
	SizeEllipsoid sizes[nSizes]; // sizes of the ellipsoids
	Overlap overlaps[nEllipsoids*(nEllipsoids - 1)/2]; // overlaps between the ellipsoids

	Par par; // parameters of the simulation

	// INITIALISATION TO NULL OF FORCES

	for (int el = 0; el < nEllipsoids; el++){
		for (int coord = 0; coord < 3; coord++){
			forces[el].el[coord] = 0;
			forces[el].dis[coord] = 0;
			forces[el].M[coord] = 0;
			forces[el].wdot[coord] = 0;
		}
	}

	// INITIALISATION OF PARAMETERS, SIZES AND ELLIPSOIDS

	init_struct(&par,&sizes[0],&ellipsoids[0]); // initialisation of the parameters, sizes and ellipsoids structures

	// INTEGRATION

	// index

	int couple; // index of the couples of particles
	int size1,size2; // index of the sizes of the particles

	// copies for modified Verlet's method

	Ellipsoid ellipsoids_[nEllipsoids]; // copies of the ellipsoids

	// measurement tools

	double dist2; // distance between the centre of the ellipsoids squared

	double p[6] = {0,0,0,0,0,0}; // momentum on each axis on 2 consecutives iterations
	double L[6] = {0,0,0,0,0,0}; // angular momentum on each axis on 2 consecutives iterations

	// initialising results file

	FILE* file = init_measure(argc,argv,&par,&sizes[0],&ellipsoids[0],&forces[0],&overlaps[0],&p[0],&L[0]); // initialises results file according to name given by user
	if (file == NULL){
        // error opening file
        printf("Error opening file for measures.");
        return 0;
    }

	// integration

	for (par.iter = 0; par.iter < par.Niter; par.iter++){

		// initialisation of the forces

		for (int el = 0; el < nEllipsoids; el++){
			for (int coord = 0; coord < 3; coord++){
				forces[el].el[coord] = 0;
				forces[el].dis[coord] = 0;
				forces[el].M[coord] = 0;
			}
		}

		// forces at t = t_0

		couple = 0; // index of the couple of particles

		for (int el1 = 0; el1 < nEllipsoids - 1; el1++){
			for (int el2 = el1 + 1; el2 < nEllipsoids; el2++){

				size1 = ellipsoids[el1].size; // index of the size of particle el1
				size2 = ellipsoids[el2].size; // index of the size of particle el2

				dist2 = 0;
				for (int coord = 0; coord < 3; coord++){
					// distance between the particles' centres squared
					dist2 += pow(ellipsoids[el1].r[coord] - ellipsoids[el2].r[coord],2);
				}

				if (dist2 < pow(sizes[size1].ax[sizes[size1].axMax] + sizes[size2].ax[sizes[size2].axMax],2)){
					// particles have a chance to overlap
					contact(&par,&overlaps[couple],&sizes[0],&ellipsoids[el1],&ellipsoids[el2]); // refreshes the overlap
				}
				else{
					// particles have no chance to overlap
					overlaps[couple].overlapping = 0;
				}

				all_forces(&par,&forces[el1],&forces[el2],&ellipsoids[el1],&ellipsoids[el2],&sizes[0],&overlaps[couple]); // forces at t = t_0

				couple++; // next couple

			}
		}

		// first increment

		for (int el = 0; el < nEllipsoids; el++){

			Euler(&ellipsoids_[el],&ellipsoids[el],&sizes[0],&forces[el]); // particle at first increment

			for (int coord = 0; coord < 3; coord++){
				// initialisation of the forces
				forces[el].el[coord] = 0;
				forces[el].dis[coord] = 0;
				forces[el].M[coord] = 0;
			}

		}

		// forces at first increment

		couple = 0; // index of the couple of particles

		for (int el1 = 0; el1 < nEllipsoids - 1; el1++){
			for (int el2 = el1 + 1; el2 < nEllipsoids; el2++){

				size1 = ellipsoids_[el1].size; // index of the size of particle el1
				size2 = ellipsoids_[el2].size; // index of the size of particle el2

				dist2 = 0;
				for (int coord = 0; coord < 3; coord++){
					// distance between the particles' centres squared
					dist2 += pow(ellipsoids_[el1].r[coord] - ellipsoids_[el2].r[coord],2);
				}

				if (dist2 < pow(sizes[size1].ax[sizes[size1].axMax] + sizes[size2].ax[sizes[size2].axMax],2)){
					contact(&par,&overlaps[couple],&sizes[0],&ellipsoids_[el1],&ellipsoids_[el2]); // refreshes the overlap
				}
				else{
					// particles have no chance to overlap
					overlaps[couple].overlapping = 0;
				}

				all_forces(&par,&forces[el1],&forces[el2],&ellipsoids_[el1],&ellipsoids_[el2],&sizes[0],&overlaps[couple]); // forces at first increment

				couple++; // next couple
			}
		}

		// modified Verlet integration

		for (int el = 0; el < nEllipsoids; el++){
			mod_Verlet(&ellipsoids[el],&ellipsoids_[el],&sizes[0],&forces[el]); // particle at t = t_0 + dt
		}

		// printing the result

		measure(file, &par, &sizes[0], &ellipsoids[0], &forces[0], &overlaps[0], &p[0], &L[0]);

	}

	fclose(file);

}
