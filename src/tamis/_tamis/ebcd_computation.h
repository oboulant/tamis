typedef struct Ebcd_Res {
   double *U;
   int n_A;
   int *A;
} Ebcd_Res;

void ebcd_compute(const double *signal, const int n_samples, const int n_dims, const double lambda, const double *weights, const double tol, Ebcd_Res *res);
