#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "maths.h"

/*pointer manipulation*/

void init(double *pointer,int k){
   // Initialises to 0 every element of list of length k pointer by pointer.
   for (int i = 0; i < k; i++){
      pointer[i] = 0;
   }
}

void equal(double *mat,double *mat_,int k){
   // Assigns mat_ to mat of length k.
   for (int i = 0; i < k; i++){
      mat[i] = mat_[i];
   }
}

/*simple mathematical operations*/

void matrix_diff(double *mat, double *a, double *b, int k){
   // Associates to mat the difference of arrays a by b of length k.
   for (int i = 0;i<k;i++){
      mat[i] = a[i] - b[i]; // difference coordinate by coordinate
   }
}

double matrix_norm(double *a, int k){
   // Returns the norm of a vector represented by array a of length k.

   double norm_ = 0; // norm of a

   for (int i = 0;i<k;i++){
      norm_ += pow(a[i],2);
   }

   norm_ = sqrt(norm_);
   return norm_;

}

void scal_prod(double *mat, double scal, double *a, int k){
   // Associates to mat the array scal * a where k is the size of a.
   for (int i = 0;i<k;i++){
      mat[i] = scal*a[i];
   }
}

void matrix_sum(double *mat, double *a, double *b, int k){
   // Assoiates to mat the sum of arrays a and b of length k.
   for (int i = 0;i<k;i++){
      mat[i] = a[i] + b[i]; // sum coordinate by coordinate
   }
}

void matrix_prod(double *mat, double *mat1, double *mat2, int i1, int j1, int j2){
   // Associates to mat the matrix product of mat1 and mat2, with i1 the number of lines of mat1 and
   // j1 and j2 the number of columns of mat1 and mat2.
   // Notice that the number of lines of mat2 have to be i1.

   for (int i = 0; i < i1; i++){
      for (int j = 0; j < j2; j++){
         mat[i * j2 + j] = 0;
         for (int k = 0; k < j1; k++){
            mat[i * j2 + j] += mat1[i * j1 + k] * mat2[k * j2 + j];
         }
      }
   }

}

void matrix_prod_333(double *mat, double *mat1, double *mat2){
   // identical to matrix_prod(mat,mat1,mat2,3,3,3)

   mat[0] = mat1[0]*mat2[0] + mat1[1]*mat2[3] + mat1[2]*mat2[6];
   mat[1] = mat1[0]*mat2[1] + mat1[1]*mat2[4] + mat1[2]*mat2[7];
   mat[2] = mat1[0]*mat2[2] + mat1[1]*mat2[5] + mat1[2]*mat2[8];

   mat[3] = mat1[3]*mat2[0] + mat1[4]*mat2[3] + mat1[5]*mat2[6];
   mat[4] = mat1[3]*mat2[1] + mat1[4]*mat2[4] + mat1[5]*mat2[7];
   mat[5] = mat1[3]*mat2[2] + mat1[4]*mat2[5] + mat1[5]*mat2[8];

   mat[6] = mat1[6]*mat2[0] + mat1[7]*mat2[3] + mat1[8]*mat2[6];
   mat[7] = mat1[6]*mat2[1] + mat1[7]*mat2[4] + mat1[8]*mat2[7];
   mat[8] = mat1[6]*mat2[2] + mat1[7]*mat2[5] + mat1[8]*mat2[8];
}

void matrix_prod_331(double *mat, double *mat1, double *mat2){
   // identical to matrix_prod(mat,mat1,mat2,3,3,1)

   mat[0] = mat1[0]*mat2[0] + mat1[1]*mat2[1] + mat1[2]*mat2[2];
   mat[1] = mat1[3]*mat2[0] + mat1[4]*mat2[1] + mat1[5]*mat2[2];
   mat[2] = mat1[6]*mat2[0] + mat1[7]*mat2[1] + mat1[8]*mat2[2];
}

void matrix_prod_313(double *mat, double *mat1, double *mat2){
   // identical to matrix_prod(mat,mat1,mat2,3,1,3)

   mat[0] = mat1[0]*mat2[0];
   mat[1] = mat1[0]*mat2[1];
   mat[2] = mat1[0]*mat2[2];

   mat[3] = mat1[1]*mat2[0];
   mat[4] = mat1[1]*mat2[1];
   mat[5] = mat1[1]*mat2[2];

   mat[6] = mat1[2]*mat2[0];
   mat[7] = mat1[2]*mat2[1];
   mat[8] = mat1[2]*mat2[2];
}

void diag_mat(double *mat, double *scal, int len){
   // Associates mat to the len*len matrix of diagonal items scal.

   for (int i = 0; i < len; i++){
      for (int j = 0; j < len; j++){
         if (i == j){
            mat[i * len + j] = scal[i];
         }
         else {
            mat[i * len + j] = 0;
         }
      }
   }
}

void transpose(double *matT, double *mat, int i, int j){
   // Associates to matT the transpose of the matrix mat with i lnes and j colums.
   for (int i_ = 0; i_ < i; i_++){
      for (int j_ = 0; j_ < j; j_++){
         matT[j_ * i + i_] = mat[i_ * j + j_];
      }
   }
}

void det2_mat(double *det, double *mat){
   // Associates det to the determinant of the 2*2 matrix mat.

   det[0] = mat[0]*mat[3] - mat[1]*mat[2]; // det = ad - bc
}

void det3_mat(double *det, double *mat){
   // Associates det to the determinant of the 3*3 matrix mat.
   // This function uses the Sarrus rule.

   det[0] = mat[0]*mat[4]*mat[8] + mat[3]*mat[7]*mat[2] + mat[6]*mat[1]*mat[5];
   det[0] -= mat[6]*mat[4]*mat[2] + mat[0]*mat[7]*mat[5] + mat[3]*mat[1]*mat[8];
}

void adj_mat3(double *adj, double *mat){
   // Associates to adj the ajoint matrix of mat of size 3*3.

   double minor[4]; // minor matrix
   int count; // count of items in the minor matrix

   for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){

         // determination of the minor matrix
         count = 0;
         for (int i_ = 0; i_ < 3; i_++){
            for (int j_ = 0; j_ < 3; j_++){
               if (i_ != i && j_ != j){
                  minor[count] = mat[i_ * 3 + j_];
                  count++;
               }
            }
         }

         // determination of the determinant
         det2_mat(&adj[j * 3 + i],&minor[0]);
         if ((i + j)%2 == 1){
            adj[j * 3 + i] *= -1;
         }

      }
   }

}

void cross(double *mat, double *vec1, double *vec2){
   // Associates to mat the cross product of R^3 vectors vec1 and vec2.
   mat[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
   mat[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
   mat[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

double angle(double *vec1, double *vec2){
   // Returns the angle between R^3 vectors vec1 and vec2.

   double cross_prod[3];
   cross(&cross_prod[0],vec1,vec2);
   double num = matrix_norm(cross_prod,3);

   double den = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];

   double ang = atan(num/den);
   return ang;

}

int max(double *list, int k){
   // Returns the index of the maximum of list of doubles of length k.

   int max_ = 0; // maximum of list
   for (int i = 1; i < k; i++){
      if (list[max_] < list[i]){
         max_ = i;
      }
   }

   return max_;

}

int min(double *list, int k){
	// Returns the index of the minimum of list of doubles of length k.

	int min_ = 0; // minimum of list
	for (int i = 1; i < k; i++){
	  if (list[min_] > list[i]){
	     min_ = i;
	  }
	}

   return min_;

}

/*operations in E^3*/

void rotation(double *Q, double *q){
   // Associates to Q the R^3 rotation matrix associated to the unit quaternion q.

   double mat0[9] = {q[3],-q[2],q[1],q[2],q[3],-q[0],-q[1],q[0],q[3]}; // \vec{q}\times. + q_0I_3 ...

   double mat1[9];  // ... squared
   matrix_prod_333(&mat1[0],&mat0[0],&mat0[0]);

   double mat2[9];
   matrix_prod_313(&mat2[0],q,q); // \vec{q}\vec{q}^T

   matrix_sum(Q,&mat1[0],&mat2[0],9); // rotation matrix in R^3

}

/*operations on quaternions*/

void qprod(double *qp, double *q1, double *q2){
   // Associates to qp the product of quaternions q1 and q2.
   qp[3] = q1[3]*q2[3] - q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2];
   qp[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
   qp[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
   qp[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
}

void qconj(double *qc, double *q){
   // Associates to qc the conjugate of quaternion q.
   qc[3] = q[3];
   for (int i = 0; i < 3; i++){
      qc[i] = -q[i];
   }
}

void scal2quat(double *quat, double scal){
	// Associates to quat the quaternion associated to the scalar scal.
	init(quat,4);

	quat[3] = scal;
}

void action_q(double *rot_vec, double *q, double *vec){
   // Associates to rot_vec the action of the unit quaternion q on a R^3 vector vec.

   double quat[4]; // quaternion associated to vec
   for (int coord = 0; coord < 3; coord++){
      quat[coord] = vec[coord];
   }
   quat[3] = 0;

   double qp[4]; // product of q and vec
   qprod(&qp[0],q,&quat[0]);

   double qc[4]; // conjugate of q
   qconj(&qc[0],q);

   double quat_[4]; // quaternion associated to the the action of q on vec
   qprod(&quat_[0],&qp[0],&qc[0]);
   equal(rot_vec,&quat_[0],3);

}
