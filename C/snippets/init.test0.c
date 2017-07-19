/*SIZES*/

int nSizes = 1; // number of sizes of ellipsoids

double m[] = {1}; // masses — as a 1xnSizes matrix
double R[] = {0.5}; // nominal radii — as a 1xnSizes matrix
double major[] = {1}; // coefficients of the major axis — as a 1xnSizes matrix
double minor[] = {1}; // coefficients of the minor axis — as a 1xnSizes matrix

// ELLIPSOIDS

int nEllipsoids = 2; // number of ellipsoids

int eSizes[] = {0,0}; // index of the sizes of the ellipsoids

double r[] = {0,0,0,2,0.5,0}; // initial positions — as a 3xnEllipsoids matrix
double v[] = {0,0,0,-0.01,0,0}; // initial velocities — as a 3xnEllipsoids matrix
double q[] = {0,0,0,1,0,0,0,1}; // initial quaternions — as a 4xnEllipsoids matrix
double w[] = {0,0,0,0,0,0}; // initial angular velocities — as a 3xnEllipsoids matrix
