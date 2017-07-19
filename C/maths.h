#ifndef MATH_H
#define MATH_H

/*PROTOTYPES*/

/*pointer manipulation*/

void init(double *pointer,int k);
   /*Initialises to 0 every element of list of length k pointer by pointer.*/

void equal(double *mat,double *mat_,int k);
   /*Assigns mat_ to mat of length k.*/

/*simple mathematical operations*/

void matrix_diff(double *mat, double *a, double *b, int k);
   /*Associates to mat the difference of arrays a by b of length k.*/

double matrix_norm(double *a, int k);
   /*Returns the norm of a vector represented by array a of length k.*/

void scal_prod(double *mat, double scal, double *a, int k);
   /*Associates to mat the array scal * a where k is the size of a.*/

void matrix_sum(double *mat, double *a, double *b, int k);
   /*Assoiates to mat the sum of arrays a and b of length k.*/

void matrix_prod(double *mat, double *mat1, double *mat2, int i1, int j1, int j2);
   /*Associates to mat the matrix product of mat1 and mat2, with i1 the number of lines of mat1 and
   j1 and j2 the number of columns of mat1 and mat2.
   Notice that the number of lines of mat2 have to be i1.*/

void matrix_prod_333(double *mat, double *mat1, double *mat2);
   /*identical to matrix_prod(mat,mat1,mat2,3,3,3)*/

void matrix_prod_331(double *mat, double *mat1, double *mat2);
   /*identical to matrix_prod(mat,mat1,mat2,3,3,1)*/

void matrix_prod_313(double *mat, double *mat1, double *mat2);
   /*identical to matrix_prod(mat,mat1,mat2,3,1,3)*/

void diag_mat(double *mat, double *scal, int len);
   /*Associates mat to the len*len matrix of diagonal items scal.*/

void transpose(double *matT, double *mat, int i, int j);
   /*Associates to matT the transpose of the matrix mat with i lnes and j colums.*/

void det2_mat(double *det, double *mat);
   /*Associates det to the determinant of the 2*2 matrix mat.*/

void det3_mat(double *det, double *mat);
   /*Associates det to the determinant of the 3*3 matrix mat.
   This function uses the Sarrus rule.*/

void adj_mat3(double *adj, double *mat);
   /*Associates to adj the ajoint matrix of mat of size 3*3.*/

void cross(double *mat, double *vec1, double *vec2);
   /*Associates to mat the cross product of R^3 vectors vec1 and vec2.*/

double angle(double *vec1, double *vec2);
   /*Returns the angle between R^3 vectors vec1 and vec2.*/

int max(double *list, int k);
   /*Returns the index of the maximum of list of doubles of length k.*/

int min(double *list, int k);
	/*Returns the index of the minimum of list of doubles of length k.*/

/*operations in E^3*/

void rotation(double *Q, double *q);
   /*Associates to Q the R^3 rotation matrix associated to the unit quaternion q.*/

/*operations on quaternions*/

void qprod(double *qp, double *q1, double *q2);
   /*Associates to qp the product of quaternions q1 and q2.*/

void qconj(double *qc, double *q);
   /*Associates to qc the conjugate of quaternion q.*/

void scal2quat(double *quat, double scal);
	/*Associates to quat the quaternion associated to the scalar scal.*/

void action_q(double *rot_vec, double *q, double *vec);
   /*Associates to rot_vec the action of the unit quaternion q on a R^3 vector vec.*/

#endif
