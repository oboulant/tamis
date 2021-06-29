cimport tamis._tamis.ebcd as ebcd

cimport cython
import numpy as np

DEFAULT_GAIN_TOL = 1e-8

cpdef ebcd_core(double[:, :] signal, double[:] weights, double lbd, double tol=DEFAULT_GAIN_TOL):
    cdef:
        int n_samples = signal.shape[0]
        int n_dims = signal.shape[1]
        ebcd.Ebcd_Res res
        int n_bkps
        int A
        cdef double[:, ::1] signal_arr = np.ascontiguousarray(signal)
        cdef double[::1] weights_arr = np.ascontiguousarray(weights)

    ebcd.ebcd_compute(&signal_arr[0, 0], n_samples, n_dims, lbd, &weights_arr[0], tol, &res)
    A_list = [el for el in res.A[:res.n_A]]
    A_arr = np.sort(A_list)
    A_arr = np.append(A_arr, [signal.shape[0]-1])
    U_list = [el for el in res.U[:n_samples*n_dims]]
    U_arr = np.array(U_list).reshape((n_samples, n_dims))
    return res.n_A, A_arr, U_arr
