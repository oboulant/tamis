cdef extern from "ebcd_computation.h":
    ctypedef struct Ebcd_Res:
        double *U
        int n_A
        int *A
    void ebcd_compute(const double *signal, const int n_samples, const int n_dims, const double lbd, const double *weights, const double tol, Ebcd_Res *res)


