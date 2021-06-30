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
n_bkps, sigma = 3, 1  # number of change points, noise standard deviation
signal, bkps = rpt.pw_constant(n, dim, n_bkps, noise_std=sigma)

#######################################
#   Compute change points
#######################################
weights = get_default_weights(n)
n_A, A, U = tam.ebcd(signal, weights, 60.0)


#######################################
#   Communicate on results
#######################################
print(f"There are {n_A} break points are :")
print(A)

rpt.display(signal, bkps, A)
plt.show()
