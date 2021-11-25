import tamis as tam
import numpy as np
import ruptures as rpt
import matplotlib.pylab as plt


def get_default_weights(n):
    i = np.arange(1, n)
    return np.sqrt(n * 1.0 / (i * 1.0 * (n - i)))


#######################################
#   Generate signal
#######################################
n, dim = 500, 3  # number of samples, dimension
n_bkps, sigma = 3, 0.2  # number of change points, noise standard deviation
signal, bkps = rpt.pw_constant(n, dim, n_bkps, noise_std=sigma, seed=1234)

#######################################
#   Compute change points
#######################################
weights = get_default_weights(n)
regularization_lambda = 30.0
n_A, A, U = tam.ebcd(signal, weights, regularization_lambda)


#######################################
#   Communicate on results
#######################################
print(f"{n_A} break points found :")
print(list(A))
print(f"True break points are :")
print(bkps)

rpt.display(signal, bkps, A)
plt.show()
