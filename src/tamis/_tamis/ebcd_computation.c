#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ebcd_computation.h"

#define DEFAULT_GAIN_TOL 1e-8
#define DEFAULT_BETA_TOL 1e-12
#define DEFAULT_MAX_IT_OUTER 1e5
#define DEFAULT_MAX_IT 1e5

/******************************************
 *
 *  Mathematical functions
 *
 ******************************************/

static inline int max_i(int num1, int num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

static inline int min_i(int num1, int num2)
{
    return (num1 > num2 ) ? num2 : num1;
}

static inline double max_f(double num1, double num2)
{
    return (num1 > num2 ) ? num1 : num2;
}

static inline double frobenius_norm(const double *Y, const int n_samples, const int n_dims)
{
    int i, j;
    double res;

    res = 0.0;

    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            res += pow(Y[i*n_dims + j], 2);
        }
    }

    return sqrt(res);
}

static inline void center_signal(const double *signal, const int n_samples, const int n_dims, double **res)
{
    int i, j;
    double *col_sum;

    col_sum = (double*)malloc(n_dims * sizeof(double));

    // Compute sum for each column
    for (j=0 ; j<n_dims ; j++)
    {
        // Initialize to 0.0
        col_sum[j] = 0.0;

        for (i=0 ; i<n_samples ; i++)
        {
            col_sum[j] += signal[i*n_dims + j];
        }
        // Divide by the number of sample to get the empirical mean
        col_sum[j] = col_sum[j] / (n_samples*1.0);
        // Substract column mean
        for (i=0 ; i<n_samples ; i++)
        {
            (*res)[i*n_dims + j] = signal[i*n_dims + j] - col_sum[j];
        }
    }

    free(col_sum);

    return;
}

static inline void cumsum(const double *signal, const int n_samples, const int n_dims, double **res)
{
    int i,j;

    // Initialize
    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] = 0.0;
        }
    }

    // Handle first line
    for (j=0 ; j<n_dims ; j++)
    {
        (*res)[j] = signal[j];
    }
    // Handle rest of the lines
    for (i=1 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] += (*res)[(i-1)*n_dims + j] + signal[i*n_dims + j];
        }
    }

    return;
}

/******************************************
 *
 *  Array logics
 *
 ******************************************/

static inline int max_index_array(const double *Y, const int n)
{
    int i, res;
    double max;

    if (n<1)
    {
        return -1;
    }
    if (n == 1)
    {
        return 0;
    }

    // Initialize
    res = 0;
    max = Y[0];

    // Scan the rest of the array
    for (i=1 ; i<n ; i++)
    {
        if (Y[i] > max)
        {
            res = i;
            max = Y[i];
        }
    }
    return res;
}

static inline void remove_element_i(const int *array_from, const int n_el, const int i_remove, int **array_to)
{
    int i, c;
    c = 0;

    for (i=0 ; i<n_el ; i++)
    {
        if (i == i_remove)
            continue;
        (*array_to)[c] = array_from[i];
        c++;
    }
}

static inline int check_i_in_A_remove_inplace(int **A, const int array_size, const int val)
{
    int i, j;

    for (i=0 ; i<array_size ; i++)
    {
        if ((*A)[i] == val)
        {
            for (j=i ; j<array_size-1 ; j++)
            {
                (*A)[j] = (*A)[j+1];
            }
            (*A)[array_size-1] = '\0';
            return 1;
        }
    }
    return 0;
}

static inline void add_element_in_A_inplace(int **A, const int max_array_size, const int array_size, const int val)
{
    if (array_size >= max_array_size)
    {
        return;
    }
    (*A)[array_size] = val;
    (*A)[array_size+1] = '\0';
    return;
}

/******************************************
 *
 *  Sparse operation involving X
 *
 ******************************************/

static inline void multiplyXnotcenteredbysparse(const double *beta, const int n_samples, const int n_dims, const double *weights, double **res)
{
    int i, j;
    double *beta_weighted;
    double *second_line;

    beta_weighted = (double*)malloc((n_samples-1) * n_dims * sizeof(double));

    // Initialize
    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] = 0.0;
        }
    }

    // Multiply each row by corresponding weights
    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            beta_weighted[i*n_dims + j] = beta[i*n_dims + j] * weights[i];
        }
    }

    // First line of the result is always filled with 0.0
    second_line = (*res)+n_dims;
    cumsum(beta_weighted, n_samples-1, n_dims, &second_line);

    free(beta_weighted);

    return;
}

void leftmultiplybyXt(const double *Y, const int n_samples, const int n_dims, const double *weights, double **res)
{
    double *Y_cumsum;
    int i, j;

    Y_cumsum = (double*)malloc(n_samples * n_dims * sizeof(double));
    cumsum(Y, n_samples, n_dims, &Y_cumsum);

    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] = weights[i] * ((i+1.0)/n_samples*Y_cumsum[(n_samples-1)*n_dims + j] - Y_cumsum[i*n_dims + j]);
        }
    }

    free(Y_cumsum);

}

void XtX(const int *A_indexes, const int A_size, const int *B_indexes, const int B_size, const double *weights, const int n_samples, double **res)
{
    int i, j;
    int u, v;

    for (i=0 ; i<A_size ; i++)
    {
        for (j=0 ; j<B_size ; j++)
        {
            u = min_i(A_indexes[i], B_indexes[j]);
            v = max_i(A_indexes[i], B_indexes[j]);
            (*res)[i*B_size + j] = weights[u] * weights[v] * (u+1.0) * (n_samples * 1.0 - (v+1.0)) / (n_samples * 1.0);
        }
    }
}

void multiplyXtXbysparse(const int *A_indexes, const int A_size, const int n_samples, const int n_dims, const double *beta, const double *weights, double **res)
{
    int i, j;
    double *beta_c, *s, *u, *pre_res;

    // If A is the empty set, results is just filled with 0.0
    if (A_size == 0)
    {
        for (i=0 ; i<n_samples-1 ; i++)
        {
            for (j=0 ; j<n_dims ; j++)
            {
                (*res)[i*n_dims + j] = 0.0;
            }
        }
        return;
    }

    // Allocate local variables
    beta_c = (double*)malloc((n_samples-1)*n_dims*sizeof(double));
    s = (double*)malloc((n_samples-1)*n_dims*sizeof(double));
    u = (double*)malloc(n_dims*sizeof(double));
    pre_res = (double*)malloc((n_samples-1)*n_dims*sizeof(double));

    // Initialize
    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            beta_c[i*n_dims + j] = 0.0;
            s[i*n_dims + j] = 0.0;
            pre_res[i*n_dims + j] = 0.0;
        }
    }
    for (j=0 ; j<n_dims ; j++)
    {
        u[j] = 0.0;
    }

    // Multiply beta by corresponding weight
    for (i=0 ; i<A_size ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            beta_c[A_indexes[i]*n_dims + j] = beta[A_indexes[i]*n_dims + j] * weights[A_indexes[i]];
        }
    }

    // Induction in reverse order of beta_c
    for (j=0 ; j<n_dims ; j++)
    {
        s[(n_samples-2)*n_dims + j] = beta_c[(n_samples-2)*n_dims + j];
    }
    for (i=1 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            s[(n_samples-2 -i)*n_dims + j] = s[(n_samples-2 - i+1)*n_dims + j] + beta_c[(n_samples-2 -i)*n_dims + j];
        }
    }

    // Compute 1 x p vector u
    for (i=0 ; i<A_size ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            u[j] += (A_indexes[i]*1.0+1.0)*beta_c[A_indexes[i]*n_dims + j];
        }
    }
    // Divide by the number of samples
    for (j=0 ; j<n_dims ; j++)
    {
        u[j] = u[j] / (n_samples*1.0);
    }

    // Handle first line of the pre-result
    for (j=0 ; j<n_dims ; j++)
    {
        pre_res[j] = s[j] - u[j];
    }
    // Handle rest of the lines
    for (i=1 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            pre_res[i*n_dims + j] = pre_res[(i-1)*n_dims + j] + s[i*n_dims + j] - u[j];
        }
    }

    // Just apply the weights
    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            (*res)[i*n_dims + j] = pre_res[i*n_dims + j] * weights[i];
        }
    }

    free(beta_c);
    free(s);
    free(u);
    free(pre_res);
}


void ebcd_compute(const double *signal, const int n_samples, const int n_dims, const double lambda, const double *weights, const double tol, Ebcd_Res *res)
{
    int i, j, n_A, p, q, n_A_not_indexes;
    int it, A_idx;
    int *A, *Ai, *A_not_indexes;
    int global_sol, it_counter, i_max_norm_A_not_indexes;
    double tol_c, lagr;
    double XitX_dot_beta_Ai, gammai;
    double *centered_signal, *beta, *C, *U;
    double *gain, *S_i, *XitX, *new_beta_i, *S, *normS;
    double temp_d, max_norm_A_not_indexes;
    double *temp_d_array;
    double *gamma;

    tol_c = tol;
    if (tol < 0.0)
    {
        tol_c = DEFAULT_GAIN_TOL;
    }

    /*
    *   Center signal
    */
    centered_signal = (double*)malloc(n_samples * n_dims * sizeof(double));
    center_signal(signal, n_samples, n_dims, &centered_signal);

    /*
    * Initialize A and beta
    */
    n_A = 0; // We start with an empty active set
    A = (int*)malloc((n_samples-1) * sizeof(int));
    A_not_indexes = (int*)malloc((n_samples-1) * sizeof(int));
    if (A != NULL)
    {
        A[0] = '\0';
    }
    beta = (double*)malloc((n_samples-1) * n_dims * sizeof(double));
    for (i=0 ; i<n_samples-1 ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            beta[i*n_dims + j] = 0.0;
        }
    }

    /*
    * Compute C = X'*Y
    */
    C = (double*)malloc((n_samples-1) * n_dims * sizeof(double));
    leftmultiplybyXt(centered_signal, n_samples, n_dims, weights, &C);

    // Stopping criteria
    global_sol = 0;
    it_counter = 0;

    while (global_sol == 0 && it_counter < DEFAULT_MAX_IT_OUTER)
    {
        // Debug
        it_counter+=1;

        // Allocate variables used in the inner for loop
        Ai = (int*)malloc((n_A-1) * sizeof(int));
        S_i = (double*)malloc(n_dims * sizeof(double));
        new_beta_i = (double*)malloc(n_dims * sizeof(double));
        temp_d_array = (double*)malloc(1 * n_dims * sizeof(double));
        XitX = (double*)malloc(1 * (n_A-1) * sizeof(double));
        gain = (double*)malloc(n_A * sizeof(double));

        // Initialize gains
        for (i=0 ; i<n_A ; i++)
        {
            gain[i] = 2.0*tol_c;
        }

        for (it=0 ; it<DEFAULT_MAX_IT ; it++)
        {
            // Update beta

            if (n_A == 0)
            {
                // If A is the empty set, then no need to update beta
                break;
            }

            A_idx = it % n_A;
            i = A[A_idx];
            remove_element_i(A, n_A, A_idx, &Ai);

            // Initialize S with C
            for (p=0 ; p<n_dims ; p++)
            {
                S_i[p] = C[i*n_dims + p];
            }
            // Then, if the remaining active set is non empty, compute corrective term of S
            if (n_A > 1)
            {
                XtX(&i, 1, Ai, n_A-1, weights, n_samples, &XitX);
                for (p=0 ; p<n_dims ; p++)
                {
                    XitX_dot_beta_Ai = 0.0;
                    for (q=0 ; q<n_A-1 ; q++)
                    {
                        XitX_dot_beta_Ai += XitX[q] * beta[Ai[q]*n_dims + p];
                    }
                    S_i[p] -= XitX_dot_beta_Ai;
                }
            }

            // Compute new beta values
            gammai = (i + 1.0) * (n_samples - i - 1.0) * pow(weights[i], 2) / (n_samples*1.0);
            temp_d = max_f(1.0 - lambda / frobenius_norm(S_i, 1, n_dims), 0.0);
            for (p=0 ; p<n_dims ; p++)
            {
                new_beta_i[p] = temp_d * S_i[p] / gammai;
                temp_d_array[p] = beta[i*n_dims + p] - new_beta_i[p];
            }
            gain[A_idx] = frobenius_norm(temp_d_array, 1, n_dims);
            for (p=0 ; p<n_dims ; p++)
            {
                beta[i*n_dims + p] = new_beta_i[p];
            }

            // Evaluate gain and stop updating beta is under a given threshold
            if (gain[max_index_array(gain, n_A)] < tol_c)
            {
                break;
            }
        }
        free(gain);
        free(XitX);
        free(temp_d_array);
        free(new_beta_i);
        free(S_i);
        free(Ai);

        // Remove from active set the zero coefficients
        for (i=(n_samples-2) ; i>=0 ; i--)
        {
            temp_d = 0.0;
            for (j=0 ; j<n_dims ; j++)
            {
                temp_d += fabs(beta[i*n_dims + j]);
            }
            if (temp_d < DEFAULT_BETA_TOL)
            {
                for (j=0 ; j<n_dims ; j++)
                {
                    beta[i*n_dims + j] = 0.0;
                }
                if (check_i_in_A_remove_inplace(&A, n_A, i) == 1)
                {
                    n_A -= 1;
                }
            }
        }

        // Check optimality

        // Allocate local variables
        temp_d_array = (double*)malloc((n_samples-1) * n_dims * sizeof(double));
        S = (double*)malloc((n_samples-1) * n_dims * sizeof(double));
        normS = (double*)malloc((n_samples-1) * sizeof(double));

        // Compute X.T @ X @ beta for corrective term of S
        multiplyXtXbysparse(A, n_A, n_samples, n_dims, beta, weights, &temp_d_array);
        // Compute S and norm of each line of S
        for (i=0 ; i<n_samples-1 ; i++)
        {
            normS[i] = 0.0;
            for (j=0 ; j<n_dims ; j++)
            {
                S[i*n_dims + j] = C[i*n_dims + j] - temp_d_array[i*n_dims + j];
                normS[i] += pow(S[i*n_dims + j ], 2);
            }
        }
        free(temp_d_array);
        free(S);

        if (n_A > 0)
        {
            // At optimality we must have normS(i)=lambda^2 for i in AS and
            // normS(i)<lambda^2 for i not in AS.
            lagr =0.0;
            for (i=0 ; i<n_A ; i++)
            {
                if (normS[A[i]] > lagr)
                {
                    lagr = normS[A[i]];
                }
            }
            if (pow(lambda, 2) < lagr)
            {
                lagr = pow(lambda, 2);
            }

            // Construct an array with all indexes not in A
            // Start with all indexes
            n_A_not_indexes = n_samples-1;
            for (i=0 ; i<n_samples-1 ; i++)
            {
                A_not_indexes[i] = i;
            }
            // Removes index
            for (i=0 ; i<n_A ; i++)
            {
                if(check_i_in_A_remove_inplace(&A_not_indexes, n_A_not_indexes, A[i]) == 1)
                {
                    // Should be the case n_A times
                    n_A_not_indexes -= 1;
                }
            }

            // Look for max of normS for indexes in A_not_indexes
            max_norm_A_not_indexes = 0.0;
            i_max_norm_A_not_indexes = -1;
            for (i=0 ; i<n_A_not_indexes ; i++)
            {
                if (normS[A_not_indexes[i]] > max_norm_A_not_indexes)
                {
                    max_norm_A_not_indexes = normS[A_not_indexes[i]];
                    i_max_norm_A_not_indexes = A_not_indexes[i];
                }
            }
            if ((n_A_not_indexes == 0) || (max_norm_A_not_indexes < lagr + tol_c))
            {
                // Optimality conditions are fulfilled: we have found the global
                // solution
                global_sol = 1;
            }
            else
            {
                // Otherwise we add the block that violates most the optimality
                // condition
                add_element_in_A_inplace(&A, n_samples-1, n_A, i_max_norm_A_not_indexes);
                n_A += 1;
                for (j=0 ; j<n_dims ; j++)
                {
                    beta[i_max_norm_A_not_indexes*n_dims + j] = 0.0;
                }
            }
        }
        else
        {
            i = max_index_array(normS, n_samples);
            if (normS[i] < pow(lambda, 2) + tol_c)
            {
                global_sol = 1;
            }
            else
            {
                add_element_in_A_inplace(&A, n_samples-1, n_A, i);
                n_A = 1;
                for (j=0 ; j<n_dims ; j++)
                {
                    beta[i*n_dims + j] = 0.0;
                }
            }
        }
        free(normS);
    }

    // Reconstruct U

    U = (double*)malloc(n_samples * n_dims * sizeof(double));
    gamma = (double*)malloc(n_dims * sizeof(double));

    // Compute X @ beta
    multiplyXnotcenteredbysparse(beta, n_samples, n_dims, weights, &U);

    // Initialize gamma
    for (j=0 ; j<n_dims ; j++)
    {
        gamma[j] = 0.0;
    }
    // Compute gamma
    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            gamma[j] += (signal[i*n_dims + j] - U[i*n_dims + j]) / (n_samples*1.0) ;
        }
    }
    // Update U with gamma
    for (i=0 ; i<n_samples ; i++)
    {
        for (j=0 ; j<n_dims ; j++)
        {
            U[i*n_dims+j] = gamma[j] + U[i*n_dims+j];
        }
    }

    // Return
    res->n_A = n_A;
    res->U = (double*)malloc(n_samples * n_dims * sizeof(double));
    res->A = (int*)malloc(n_A * sizeof(int));
    memcpy(res->U, U, n_samples*n_dims*sizeof(double));
    memcpy(res->A, A, n_A*sizeof(int));



    free(gamma);
    free(U);
    free(C);
    free(beta);
    free(A_not_indexes);
    free(A);
    free(centered_signal);

    return;
}
