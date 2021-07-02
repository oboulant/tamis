#include <stdlib.h>
#include <math.h>

#include "unity.h"
#include "ebcd_computation.c" // Include directly .c because of static inline

void setUp(void) {
    // set stuff up here
}

void tearDown(void) {
    // clean stuff up here
}

void test_cumsum(void)
{

    static const double signal[6]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double *res;

    static double expected1[6]={1.0, 2.0, 3.0, 5.0, 7.0, 9.0};
    static double expected2[6]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    static double expected3[6]={1.0, 3.0, 6.0, 10.0, 15.0, 21.0};

    res = (double*)malloc(6 * sizeof(double));

    cumsum(signal, 2, 3, &res);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected1, res, 6);

    cumsum(signal, 1, 6, &res);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected2, res, 6);

    cumsum(signal, 6, 1, &res);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected3, res, 6);

    free(res);

}

void test_frobenius_norm(void)
{
    static const double signal[6]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    double res;

    res = frobenius_norm(signal, 2, 3);
    TEST_ASSERT_EQUAL_DOUBLE(sqrt(91.0), res);

    res = frobenius_norm(signal, 1, 6);
    TEST_ASSERT_EQUAL_DOUBLE(sqrt(91.0), res);

    res = frobenius_norm(signal, 6, 1);
    TEST_ASSERT_EQUAL_DOUBLE(sqrt(91.0), res);

}

void test_center_signal(void)
{
    static const double signal[6]={1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    double *res;

    static double expected1[6]={-1.5, -1.5, -1.5, 1.5, 1.5, 1.5};
    static double expected2[6]={-2.0, -2.0, 0.0, 0.0, 2.0, 2.0};
    static double expected3[6]={-2.5, -1.5, -0.5, 0.5, 1.5, 2.5};
    static double expected4[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    res = (double*)malloc(6 * sizeof(double));

    center_signal(signal, 2, 3, &res);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected1, res, 6);

    center_signal(signal, 3, 2, &res);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected2, res, 6);

    center_signal(signal, 6, 1, &res);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected3, res, 6);

    center_signal(signal, 1, 6, &res);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected4, res, 6);

    free(res);

}

void test_max_index_array(void)
{
    static const double signal[6]={1.0, 2.0, 0.5, 4.0, 12.0, 6.0};
    static const double signal2[6]={12.0, 2.0, 0.5, 4.0, 12.0, 6.0};
    static const double signal3[6]={1.0, 2.0, 0.5, 4.0, 5.0, 6.0};

    TEST_ASSERT_EQUAL_INT(4, max_index_array(signal, 6));
    TEST_ASSERT_EQUAL_INT(0, max_index_array(signal2, 6));
    TEST_ASSERT_EQUAL_INT(5, max_index_array(signal3, 6));

    TEST_ASSERT_EQUAL_INT(-1, max_index_array(signal, -1));

    TEST_ASSERT_EQUAL_INT(0, max_index_array(signal, 1));

}

void test_remove_element_i(void)
{
    static const int array_from[6]={2, 7, 2, 9, 12, 6};

    static const int expected1[5]={2, 7, 2, 12, 6};
    static const int expected2[5]={7, 2, 9, 12, 6};
    static const int expected3[5]={2, 7, 2, 9, 12};

    int *res;

    res = (int*)malloc(5 * sizeof(int));

    remove_element_i(array_from, 6, 3, &res);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected1, res, 5);

    remove_element_i(array_from, 6, 0, &res);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected2, res, 5);

    remove_element_i(array_from, 6, 5, &res);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected3, res, 5);

    remove_element_i(array_from, 6, 16, &res);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected3, res, 5);

    remove_element_i(array_from, 6, -2, &res);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected3, res, 5);

    free(res);
}

void get_array_example(int ** array)
{
    (*array)[0]=2;
    (*array)[1]=7;
    (*array)[2]=2;
    (*array)[3]=9;
    (*array)[4]=12;
    (*array)[5]=6;
}

void test_check_i_in_A_remove_inplace(void)
{
    int *array;
    int res;
    static const int expected1[5]={7, 2, 9, 12, 6};
    static const int expected2[6]={2, 7, 2, 9, 12, 6};
    static const int expected3[5]={2, 7, 2, 9, 12};

    array = (int*)malloc(6 * sizeof(int));

    get_array_example(&array);
    res = check_i_in_A_remove_inplace(&array, 6, 2);
    TEST_ASSERT_EQUAL_INT(1, res);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected1, array, 5);
    TEST_ASSERT_EQUAL_CHAR ('\0', array[5]);

    get_array_example(&array);
    res = check_i_in_A_remove_inplace(&array, 6, 24);
    TEST_ASSERT_EQUAL_INT(0, res);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected2, array, 6);

    get_array_example(&array);
    res = check_i_in_A_remove_inplace(&array, 6, 6);
    TEST_ASSERT_EQUAL_INT(1, res);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected3, array, 5);
    TEST_ASSERT_EQUAL_CHAR ('\0', array[5]);

    free(array);
}

void test_add_element_in_A_inplace(void)
{
    int *array;
    static const int expected1[7]={2, 7, 2, 9, 12, 6, 24};

    array = (int*)malloc(7 * sizeof(int));

    get_array_example(&array);
    add_element_in_A_inplace(&array, 7, 6, 24);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected1, array, 7);
    TEST_ASSERT_EQUAL_CHAR ('\0', array[7]);

    get_array_example(&array);
    add_element_in_A_inplace(&array, 6, 6, 24);
    TEST_ASSERT_EQUAL_INT_ARRAY(expected1, array, 6);

    free(array);
}

void get_weights(const int n, double **weights)
{
    int i;

    for (i=1 ; i<n ; i++)
    {
        (*weights)[i-1] = sqrt(1.0*n / (i*1.0*(n-i)));
    }

    return;
}

void get_X(const int n, const double *weights, double **X)
{
    int i, j;

    for (i=0 ; i<n ; i++)
    {
        for (j=0 ; j<n-1 ; j++)
        {
            if (i>j)
            {
                (*X)[i*(n-1) + j] = weights[j];
            }
            else
            {
                (*X)[i*(n-1) + j] = 0.0;
            }
        }
    }
    return;
}

void get_transpose(const double *from, const int nrows, const int ndims, double **to)
{
    int i, j;

    for (i=0 ; i<nrows ; i++)
    {
        for (j=0 ; j<ndims ; j++)
        {
            (*to)[j*nrows+i] = from[i*ndims + j];
        }
    }
    return;
}

// bkps = [12, 24, 37, 50]
const double signal_[150]={1.15790004,  -8.60728516,   6.87092146, \
        -2.80559417, -18.47842045,   1.62604215, \
        -5.49637749,  -6.37306161,  -1.51390604, \
        -9.47648535,  -8.07709835,   5.96563165, \
        -6.16348482, -12.01650786,  -0.99732217, \
        -3.14815796,  -7.38423577,   2.31966775, \
        -0.57026315, -11.80957253,   1.62585868, \
         0.50563261,  -4.00268522,   8.37553539, \
        -6.82554358,  -3.22178915,  -6.32605359, \
        -7.05602186,  -4.14762744,   4.20614725, \
        -0.60247347,  -3.10067907,  -0.88843904, \
        -7.69134562, -14.95630247,   1.78342254, \
       -14.39716714, -16.21243619,   0.76256198, \
        -9.09665492, -14.55606032,   0.43450941, \
        -6.77445738, -13.89650098,   2.69416362, \
        -5.84802438, -15.60304888,   3.58874691, \
        -3.1862639 , -14.95953651,   1.39064533, \
        -3.49322869,  -8.12888333,   8.6755258 , \
        -5.00610623,  -6.06001285,   9.04642104, \
        -2.70406868, -14.49574451,   2.99946756, \
         4.18128135, -16.68998497,   7.14525781, \
        -7.09945432, -11.97200431,   3.48567809, \
        -1.7132491 ,  -8.76909664,   6.81425253, \
        -5.57546619,  -9.00198971,   2.14694413, \
       -10.97381718, -10.62343508,   4.12042863, \
        -1.24801229,  -6.02222576,   5.0142595 , \
       -10.42589972, -15.47613235,  18.61689824, \
         6.35663608, -13.5516717 ,   2.34613854, \
        -7.13452075,  -6.12135608,   7.06133773, \
       -14.40797761,  -4.37200851,  14.41162984, \
        -7.39150824,  -8.33482232,  13.35155701, \
        -8.60764431,  -3.07581032,  -0.89752421, \
        -7.61682916,   2.79880012,   9.61382195, \
        -6.59777368, -12.8260744 ,  20.10302447, \
        -7.65154411,  -2.23356291,  12.82380499, \
        -8.42106658, -13.31546695,  10.20059699, \
        -2.71058265,   0.22698968,   3.59227641, \
        -1.52109416, -11.08335871,  -0.8474314 , \
        -9.94653212, -14.97327242,  10.08942715, \
        -4.60681859, -19.25782157,   0.10799169, \
       -12.37665121,  -9.64275777,  -1.63635742, \
       -11.09669429,  -9.72381691,   2.94584936, \
       -16.30205612,  -9.60762918,   6.27885515, \
       -12.23438873, -17.69412005,   6.04352236, \
         0.24462888, -18.45958861,   7.65107037, \
       -13.87844837, -16.97797008,   5.11657703, \
        -5.99435023, -14.70807776,  -8.07613473, \
       -12.91214113, -16.14605663,  -2.33663693, \
        -1.63708412, -15.95583197,   3.24735593, \
       -10.85954272, -23.27198132,  11.65095138};

void inner_leftmultiplybyXt(const int nsamples, const int ndims, double **res, double **expected)
{
    int i, j, k;
    double *weights, *X, *centered_X, *centered_X_T, *centered_signal;

    weights = (double*)malloc(nsamples * sizeof(double));
    X = (double*)malloc(nsamples*(nsamples-1) * sizeof(double));
    centered_X = (double*)malloc(nsamples*(nsamples-1) * sizeof(double));
    centered_X_T = (double*)malloc((nsamples-1)*nsamples * sizeof(double));
    centered_signal = (double*)malloc(nsamples * ndims * sizeof(double));

    get_weights(nsamples, &weights);
    get_X(nsamples, weights, &X);
    center_signal(X, nsamples, nsamples-1, &centered_X);
    get_transpose(centered_X, nsamples, nsamples-1, &centered_X_T);
    center_signal(signal_, nsamples, ndims, &centered_signal);
    leftmultiplybyXt(centered_signal, nsamples, ndims, weights, res);
    for (i=0 ; i<nsamples-1 ; i++)
    {
        for (j=0 ; j<ndims ; j++)
        {
            (*expected)[i*ndims + j] = 0.0;
            for (k=0 ; k<nsamples ; k++)
            {
                (*expected)[i*ndims + j] += centered_X_T[i*nsamples + k] * centered_signal[k*ndims+j];
            }
        }
    }

    free(weights);
    free(X);
    free(centered_X);
    free(centered_X_T);
    free(centered_signal);


}

void test_leftmultiplybyXt(void)
{
    double *res, *expected;
    int nsamples, ndims;

    /*
    *   50 x 3
    */
    nsamples = 50;
    ndims = 3;
    res = (double*)malloc((nsamples-1) * ndims * sizeof(double));
    expected = (double*)malloc((nsamples-1) * ndims * sizeof(double));
    inner_leftmultiplybyXt(nsamples, ndims, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, (nsamples-1)*ndims);
    free(res);
    free(expected);

    /*
    *   150 x 1
    */
    nsamples = 150;
    ndims = 1;
    res = (double*)malloc((nsamples-1) * ndims * sizeof(double));
    expected = (double*)malloc((nsamples-1) * ndims * sizeof(double));
    inner_leftmultiplybyXt(nsamples, ndims, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, (nsamples-1)*ndims);
    free(res);
    free(expected);

    /*
    *   75 x 2
    */
    nsamples = 75;
    ndims = 2;
    res = (double*)malloc((nsamples-1) * ndims * sizeof(double));
    expected = (double*)malloc((nsamples-1) * ndims * sizeof(double));
    inner_leftmultiplybyXt(nsamples, ndims, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, (nsamples-1)*ndims);
    free(res);
    free(expected);

}

void inner_XtX(const int nsamples, const int *A, const int A_size, const int *B, const int B_size, double **res, double **expected)
{
    int i, j, k;
    double *weights, *X, *centered_X;

    weights = (double*)malloc((nsamples-1) * sizeof(double));
    get_weights(nsamples, &weights);
    XtX(A, A_size, B, B_size, weights, nsamples, res);

    X = (double*)malloc(nsamples*(nsamples-1) * sizeof(double));
    centered_X = (double*)malloc(nsamples*(nsamples-1) * sizeof(double));
    get_X(nsamples, weights, &X);
    center_signal(X, nsamples, nsamples-1, &centered_X);
    for (i=0 ; i<A_size ; i++)
    {
        for (j=0 ; j<B_size ; j++)
        {
            (*expected)[i*B_size + j] = 0.0;
            for (k=0 ; k<nsamples ; k++)
            {
                (*expected)[i*B_size + j] += centered_X[k*(nsamples-1) + A[i]] * centered_X[k*(nsamples-1) + B[j]];
            }
        }
    }

    free(weights);
    free(X);
    free(centered_X);
}

void test_XtX(void)
{
    double *res, *expected;
    int A_size, B_size, nsamples;
    int *A, *B;

    /*
    *   nsamples = 50
    *   A_size = 5
    *   B_size = 7
    */
    nsamples = 50;
    A_size = 5;
    B_size = 7;
    A = (int*)malloc(A_size * sizeof(int));
    A[0] = 2;
    A[1] = 6;
    A[2] = 17;
    A[3] = 12;
    A[4] = 46;
    B = (int*)malloc(B_size * sizeof(int));
    B[0] = 8;
    B[1] = 36;
    B[2] = 27;
    B[3] = 23;
    B[4] = 48;
    B[5] = 16;
    B[6] = 3;
    res = (double*)malloc(A_size * B_size * sizeof(double));
    expected = (double*)malloc(A_size * B_size * sizeof(double));
    inner_XtX(nsamples, A, A_size, B, B_size, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, A_size*B_size);
    free(res);
    free(expected);
    free(A);
    free(B);

    /*
    *   nsamples = 50
    *   A_size = 1
    *   B_size = 7
    */
    nsamples = 50;
    A_size = 1;
    B_size = 7;
    A = (int*)malloc(A_size * sizeof(int));
    A[0] = 37;
    B = (int*)malloc(B_size * sizeof(int));
    B[0] = 8;
    B[1] = 36;
    B[2] = 27;
    B[3] = 23;
    B[4] = 48;
    B[5] = 16;
    B[6] = 3;
    res = (double*)malloc(A_size * B_size * sizeof(double));
    expected = (double*)malloc(A_size * B_size * sizeof(double));
    inner_XtX(nsamples, A, A_size, B, B_size, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, A_size*B_size);
    free(res);
    free(expected);
    free(A);
    free(B);

    /*
    *   nsamples = 50
    *   A_size = 1
    *   B_size = 1
    */
    nsamples = 50;
    A_size = 1;
    B_size = 1;
    A = (int*)malloc(A_size * sizeof(int));
    A[0] = 37;
    B = (int*)malloc(B_size * sizeof(int));
    B[0] = 22;
    res = (double*)malloc(A_size * B_size * sizeof(double));
    expected = (double*)malloc(A_size * B_size * sizeof(double));
    inner_XtX(nsamples, A, A_size, B, B_size, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, A_size*B_size);
    free(res);
    free(expected);
    free(A);
    free(B);
}

const double beta_[150]={0.0,  0.0,   0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.7, -0.4,  3.9, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.1, 8.3,  -12.4, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0, \
        0.0, 0.0,  0.0};

void inner_multiplyXtXbysparse(const int nsamples, const int ndims, const int *A, const int A_size, const double *beta, double **res, double **expected)
{
    double *weights, *X, *centered_X, *XtX;
    int i, j, k;

    weights = (double*)malloc((nsamples-1) * sizeof(double));
    get_weights(nsamples, &weights);
    multiplyXtXbysparse(A, A_size, nsamples, ndims, beta_, weights, res);

    X = (double*)malloc(nsamples*(nsamples-1) * sizeof(double));
    centered_X = (double*)malloc(nsamples*(nsamples-1) * sizeof(double));
    XtX = (double*)malloc((nsamples-1)*(nsamples-1) * sizeof(double));
    get_X(nsamples, weights, &X);
    center_signal(X, nsamples, nsamples-1, &centered_X);

    for (i=0 ; i<nsamples-1 ; i++)
    {
        for (j=0 ; j<nsamples-1 ; j++)
        {
            XtX[i*(nsamples-1) + j] = 0.0;
            for (k=0 ; k<nsamples ; k++)
            {
                XtX[i*(nsamples-1) + j] += centered_X[k*(nsamples-1) +i] * centered_X[k*(nsamples-1) +j];
            }
        }
    }
    for (i=0 ; i<nsamples-1 ; i++)
    {
        for (j=0 ; j<ndims ; j++)
        {
            (*expected)[i*ndims + j] = 0.0;
            for (k=0 ; k<nsamples-1 ; k++)
            {
                (*expected)[i*ndims + j] += XtX[i*(nsamples-1) + k] * beta[k*ndims + j];
            }
        }
    }
    free(weights);
    free(X);
    free(centered_X);
    free(XtX);
    return;
}

void test_multiplyXtXbysparse(void)
{
    int A_size, nsamples, ndims, i;
    int *A;
    double *res, *expected, *beta;

    /*
    *   nsamples = 50
    *   A_size = 2
    */
    nsamples = 50;
    ndims = 3;
    A_size = 2;
    A = (int*)malloc(A_size * sizeof(int));
    A[0] = 12;
    A[1] = 27;
    res = (double*)malloc((nsamples-1)*ndims * sizeof(double));
    expected = (double*)malloc((nsamples-1)*ndims * sizeof(double));
    beta = (double*)malloc((nsamples-1)*ndims * sizeof(double));
    memcpy(beta, beta_, (nsamples-1)*ndims*sizeof(double));
    inner_multiplyXtXbysparse(nsamples, ndims, A, A_size, beta, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, (nsamples-1)*ndims);
    free(res);
    free(expected);
    free(beta);
    free(A);

    /*
    *   nsamples = 50
    *   A_size = 0
    */
    nsamples = 50;
    ndims = 3;
    A_size = 0;
    A = (int*)malloc(A_size * sizeof(int));
    res = (double*)malloc((nsamples-1)*ndims * sizeof(double));
    expected = (double*)malloc((nsamples-1)*ndims * sizeof(double));
    beta = (double*)malloc((nsamples-1)*ndims * sizeof(double));
    memcpy(beta, beta_, (nsamples-1)*ndims*sizeof(double));
    for (i=0 ; i<ndims ; i++)
    {
        beta[12*ndims + i] = 0.0;
        beta[27*ndims + i] = 0.0;
    }
    inner_multiplyXtXbysparse(nsamples, ndims, A, A_size, beta, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, (nsamples-1)*ndims);
    free(res);
    free(expected);
    free(beta);
    free(A);

    return;
}

void inner_multiplyXnotcenteredbysparse(const double *beta, const int nsamples, const int ndims, double **res, double **expected)
{
    double *weights, *X;
    int i, j, k;

    weights = (double*)malloc((nsamples-1) * sizeof(double));
    get_weights(nsamples, &weights);
    multiplyXnotcenteredbysparse(beta, nsamples, ndims, weights, res);

    X = (double*)malloc(nsamples*(nsamples-1) * sizeof(double));
    get_X(nsamples, weights, &X);
    for (i=0 ; i<nsamples ; i++)
    {
        for (j=0 ; j<ndims ; j++)
        {
            (*expected)[i*ndims + j] = 0.0;
            for (k=0 ; k<nsamples-1 ; k++)
            {
                (*expected)[i*ndims + j] += X[i*(nsamples-1) + k] * beta[k*ndims + j];
            }
        }
    }

    free(weights);
    free(X);
}

void test_multiplyXnotcenteredbysparse(void)
{
    double *res, *expected, *beta;
    int nsamples, ndims, i;

    /*
    *   nsamples = 50
    *   ndims = 3
    *   beta
    */
    nsamples = 50;
    ndims = 3;
    beta = (double*)malloc((nsamples-1)*ndims * sizeof(double));
    memcpy(beta, beta_, (nsamples-1)*ndims*sizeof(double));
    res = (double*)malloc(nsamples*ndims * sizeof(double));
    expected = (double*)malloc(nsamples*ndims * sizeof(double));
    inner_multiplyXnotcenteredbysparse(beta, nsamples, ndims, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, nsamples*ndims);
    free(beta);
    free(res);
    free(expected);

    /*
    *   nsamples = 50
    *   ndims = 3
    *   beta all to 0.0
    */
    nsamples = 50;
    ndims = 3;
    beta = (double*)malloc((nsamples-1)*ndims * sizeof(double));
    memcpy(beta, beta_, (nsamples-1)*ndims*sizeof(double));
    for (i=0 ; i<ndims ; i++)
    {
        beta[12*ndims + i] = 0.0;
        beta[27*ndims + i] = 0.0;
    }
    res = (double*)malloc(nsamples*ndims * sizeof(double));
    expected = (double*)malloc(nsamples*ndims * sizeof(double));
    inner_multiplyXnotcenteredbysparse(beta, nsamples, ndims, &res, &expected);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(expected, res, nsamples*ndims);
    free(beta);
    free(res);
    free(expected);

    return;

}


int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_cumsum);
    RUN_TEST(test_frobenius_norm);
    RUN_TEST(test_center_signal);
    RUN_TEST(test_max_index_array);
    RUN_TEST(test_remove_element_i);
    RUN_TEST(test_check_i_in_A_remove_inplace);
    RUN_TEST(test_add_element_in_A_inplace);
    RUN_TEST(test_leftmultiplybyXt);
    RUN_TEST(test_XtX);
    RUN_TEST(test_multiplyXtXbysparse);
    RUN_TEST(test_multiplyXnotcenteredbysparse);
    return UNITY_END();
}
