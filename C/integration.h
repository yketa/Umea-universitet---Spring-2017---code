#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "param.h"
#include "ellipsoids.h"
#include "maths.h"
#include "init.h"

/*PROTOTYPES*/

void all_forces(Par *par, Forces *forces1, Forces *forces2, Ellipsoid *ellipsoid1,
	Ellipsoid *ellipsoid2, SizeEllipsoid *sizes, Overlap *overlap);
	/*Adds to forces_i the contribution of the forces and moments exerted by ellipsoid_j on ellispoid_i.
	Thanks to Newton (III), only the axis of one ellipsoid is necessary to this computation.*/

void Euler(Ellipsoid *ellipsoid_, Ellipsoid *ellipsoid, SizeEllipsoid *sizes, Forces *forces);
	/*Associates to ellipsoid_ the integration of ellipsoid to the following time step (t + dt).*/

void mod_Verlet(Ellipsoid *ellipsoid, Ellipsoid *ellipsoid_, SizeEllipsoid *sizes, Forces *forces);
	/* Associates to ellipsoid its integration to t+dt with ellipsoid_ its integration to t+dt with the
	Euler method and forces the forces and moments exerted on the ellipsoid according to ellipsoid_.*/

#endif
