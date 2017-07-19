void translation(double *T, double *u){
   // Associates to T the R^3 translation matrix associated to u.

   init(T,16);
   for (int i = 0; i < 3; i++){
      T[i*4 + 3] = u[i];
      T[i*4 + i] = 1;
   }
   T[15] = 1;

}

void dilatation(double *X, double a){
   // Associates to X the R^3 dilatation matrix of factor a.

   init(X,16);
   for (int i = 0; i < 3; i++){
      X[i*4 + i] = a;
   }
   X[15] = 1;
   
}
