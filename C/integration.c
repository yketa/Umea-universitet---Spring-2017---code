#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "integration.h"

void all_forces(Par *par, Forces *forces1, Forces *forces2, Ellipsoid *ellipsoid1,
	Ellipsoid *ellipsoid2, SizeEllipsoid *sizes, Overlap *overlap){
	// Adds to forces_i the contribution of the forces and moments exerted by ellipsoid_j on ellispoid_i.
	// Thanks to Newton (III), only the axis of one ellipsoid is necessary to this computation.

	if (!overlap->overlapping){
	  // ellipsoids are not overlapping
	  return;
	}

	// Elastic force

	double n1[3]; // non-unit surface vector of ellipsoid1 at overlap->r_c
    double diff[3]; // ellipsoid2->r - ellipsoid1->r
    double scal = 0; // diff . n1
    double fel12[3]; // elastic force exerted by ellipsoid2 on ellipsoid1

    n_surface_vec(&n1[0],ellipsoid1,sizes,overlap->r_c);

    matrix_diff(&diff[0],ellipsoid2->r,ellipsoid1->r,3);

    for (int coord = 0; coord < 3; coord++){
    	scal += diff[coord] * n1[coord];
    }

    scal_prod(&fel12[0],- par->ke * overlap->mu * (1 - overlap->mu) / scal,&n1[0],3);

	matrix_sum(forces1->el,forces1->el,&fel12[0],3);
	matrix_diff(forces2->el,forces2->el,&fel12[0],3); // fel21 = - fel12 (Newton (III))

	// Dissipative force

	double ww1[3]; // rotation vector of ellipsoid1 in the world frame
	double ww2[3]; // rotation vector of ellispoid2 in the world frame
	for (int i = 0; i < 3; i++){
		ww1[i] = 0;
		ww2[i] = 0;
		for (int j = 0; j < 3; j++){
			ww1[i] += ellipsoid1->Q[i*3 + j]*ellipsoid1->w[j]; // ww1 = ellipsoid1->Q wb1
			ww2[i] += ellipsoid2->Q[i*3 + j]*ellipsoid2->w[j]; // ww2 = ellipsoid2->Q wb2
		}
	}

	double radius1[3]; // overlap->r_c - ellipsoid1->r
	double radius2[3]; // overlap->r_c - ellipsoid2->r
	matrix_diff(&radius1[0],overlap->r_c,ellipsoid1->r,3);
	matrix_diff(&radius2[0],overlap->r_c,ellipsoid2->r,3);

	double rotv1[3]; // velocity of r_c1 in ellispoid1 frame
	double rotv2[3]; // velocity of r_c2 in ellispoid2 frame
	cross(&rotv1[0],&ww1[0],&radius1[0]);
	cross(&rotv2[0],&ww2[0],&radius2[0]);

	matrix_sum(overlap->v1c2,ellipsoid1->v,&rotv1[0],3); // local velocity of ellipsoid1 at its point of contact with ellispoid2
	matrix_sum(overlap->v2c1,ellipsoid2->v,&rotv2[0],3); // local velocity of ellipsoid2 at its point of contact with ellispoid1

	double fdis12[3]; // dissipative force exerted by ellipsoid2 on ellipsoid1
	for (int coord = 0; coord < 3; coord++){
		fdis12[coord] = - par->kd * (overlap->v1c2[coord] - overlap->v2c1[coord]);
	}
	matrix_sum(forces1->dis,forces1->dis,&fdis12[0],3);
	matrix_diff(forces2->dis,forces2->dis,&fdis12[0],3); // fdis21 = - fdis12 (Newton (III))

	// Total force

	double ftot12[3]; // total force exerted on ellipsoid1 by ellipsoid2
	matrix_sum(&ftot12[0],&fel12[0],&fdis12[0],3);

	// Moments of forces

	double M1w[3]; // moments of forces exerted on ellipsoid1 by ellipsoid2 in the world frame
	double M2w[3]; // moments of forces exerted on ellipsoid2 by ellipsoid1 in the world frame
	cross(&M1w[0],&radius1[0],&ftot12[0]); // M1 = r1 x ftot1
	cross(&M2w[0],&ftot12[0],&radius2[0]); // M2 = r2 x ftot2 = r2 x (- ftot1) = ftot1 x r2

	double M1b[3]; // moments of forces exerted on ellipsoid1 by ellipsoid2 in the body frame
	double M2b[3]; // moments of forces exerted on ellipsoid2 by ellipsoid1 in the body frame
	for (int i = 0; i < 3; i++){
		M1b[i] = 0;
		M2b[i] = 0;
		for (int j = 0; j < 3; j++){
			M1b[i] += ellipsoid1->Q[j*3 + i]*M1w[j]; // M1b = Q^-1 M1w = tQ M1w
			M2b[i] += ellipsoid2->Q[j*3 + i]*M2w[j]; // M1b = Q^-1 M1w = tQ M1w
		}
	}

	matrix_sum(forces1->M,forces1->M,&M1b[0],3);
	matrix_sum(forces2->M,forces2->M,&M2b[0],3);

}

void Euler(Ellipsoid *ellipsoid_, Ellipsoid *ellipsoid, SizeEllipsoid *sizes, Forces *forces){
	// Associates to ellipsoid_ the integration of ellipsoid to the following time step (t + dt).

	wdot(ellipsoid,sizes,forces); // rotation vector of the particle
	
	double qp[4]; // first derivative of the quaternion
	qdot(&qp[0],ellipsoid);
	double qpnorm = matrix_norm(&qp[0],4); // norm of the first derivative of the quaternion
	double corr = 0; // correction for the integration of the quaternion
	if (qpnorm != 0){
		// the quaternion is non-constant
		corr = tan(dt*qpnorm)/qpnorm;
	}
	double qnorm2 = 0; // square of the norm of the integrated quaternion

	for (int coord = 0; coord < 3; coord++){

		ellipsoid_->r[coord] = ellipsoid->r[coord] + dt*ellipsoid->v[coord]; // new position

		ellipsoid_->v[coord] = ellipsoid->v[coord] + dt*(forces->el[coord] + forces->dis[coord]); // new velocity

		ellipsoid_->w[coord] = ellipsoid->w[coord] + dt*forces->wdot[coord]; // new rotation vector

		ellipsoid_->q[coord] = ellipsoid->q[coord] + corr*qp[coord]; // new quaternion
		qnorm2 += pow(ellipsoid_->q[coord],2);

	}
	ellipsoid_->q[3] = ellipsoid->q[3] + corr*qp[3];
	qnorm2 += pow(ellipsoid_->q[3],2);
	scal_prod(ellipsoid_->q,1/sqrt(qnorm2),ellipsoid_->q,4); // renormalisation of the integrated quaternion

	ellipsoid_->size = ellipsoid->size; // index of the size of the ellipsoid

	rotation_matrix(ellipsoid_); // rotation matrix of the ellipsoid

}

void mod_Verlet(Ellipsoid *ellipsoid, Ellipsoid *ellipsoid_, SizeEllipsoid *sizes, Forces *forces){
	// Associates to ellipsoid its integration to t+dt with ellipsoid_ its integration to t+dt with the
	// Euler method and forces the forces and moments exerted on the ellipsoid according to ellipsoid_.

	double vcoord; // velocity of the particle at t = t_0 along one axis

	double wp0[3]; // first derivative of the rotation of the particle at t = t_0
	equal(&wp0[0],forces->wdot,3);
	wdot(ellipsoid_,sizes,forces); // rotation vector of the particle at first increment

	double qp0[4]; // first derivative of the quaternion of the particle at t = t_0
	qdot(&qp0[0],ellipsoid);
	double qp1[4]; // first derivative of the quaternion of the particle at first increment
	qdot(&qp1[0],ellipsoid_);
	double dq[4]; // increment of the quaternion
	double dqnorm = 0; // norm of the increment of the quaternion
	for (int coord = 0; coord < 4; coord++){
		dq[coord] = 0.5*dt*(qp0[coord] + qp1[coord]);
		dqnorm += pow(dq[coord],2);
	}
	double corr = 0; // corection for the integration of the quaternion
	if (dqnorm != 0){
		dqnorm = sqrt(dqnorm);
		corr = tan(dqnorm)/dqnorm; 
	}
	double qnorm2 = 0; // square of the norm of the integrated quaternion

	for (int coord = 0; coord < 3; coord ++){

		vcoord = ellipsoid->v[coord]; // velocity of the particle at t = t_0
		ellipsoid->v[coord] = 0.5*(ellipsoid_->v[coord] + ellipsoid->v[coord] + (dt/(sizes + ellipsoid->size)->m) * (forces->el[coord] + forces->dis[coord])); // velocity of the particle at t = t_0 + dt

		ellipsoid->r[coord] = ellipsoid->r[coord] + 0.5*dt*(vcoord + ellipsoid->v[coord]); // position of the particle at t = t_0 + dt

		ellipsoid->w[coord] = ellipsoid->w[coord] + 0.5*dt*(wp0[coord] + forces->wdot[coord]); // rotation vector of the particle at t = t_0 + dt

		ellipsoid->q[coord] = ellipsoid->q[coord] + corr*dq[coord]; // quaternion of the particle at t = t_0 + dt
		qnorm2 += pow(ellipsoid->q[coord],2);

	}
	ellipsoid->q[3] = ellipsoid->q[3] + corr*dq[3];
	qnorm2 += pow(ellipsoid->q[3],2);
	scal_prod(ellipsoid->q,1/sqrt(qnorm2),ellipsoid->q,4); // renormalisation of the integrated quaternion

	rotation_matrix(ellipsoid); // rotation matrix of the ellipsoid

}
