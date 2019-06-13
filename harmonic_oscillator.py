import numpy as np
import scipy
import matplotlib.pyplot as plt


def histogram(A, x):
    """
    p(A) is the probability distribution of x when
    x(t) = A cos(wt) , with w = sqrt(k/m)
    solution of the harmonic oscillator with amplitude A
    A is proportional to the total energy.
    E = 1/2 k A^2
    A = sqrt(2E/k)
    """
    return 1 / (np.pi * np.sqrt(1 - (x / A)**2))


def main()


k = 0.5
m = 2
kb = 1.38064852E-23  # m2 kg s-2 K-1
temp = 351
hbar = 1.0545718E-34  # kg*m^2/s
x_val = np.arange(-4, 4, 0.1)
y_val = np.zeros(len(x_val))
 f = np.zeros(len(x_val))

 # Here a sum is done over various values of A, weighted by their relative Boltzmann factor
 # Calculating partition function:
  Z = kb * temp * np.sqrt(m / k) / hbar
   for A in np.arange(0.01, 10, 0.01):
        for i, x_tmp in enumerate(x_val):
            f[i] = histogram(A, x_tmp)
            if x_tmp < 0.99 * A and x_tmp > -0.99 * A:
                # this constraint is because at -1 and 1 the function diverges
                exponent = -(0.5 * k * A**2)  # /(kb*temp)
                y_val[i] += np.exp(exponent) * f[i]
    y_val = y_val / Z
    plt.plot(x_val, y_val)
    gauss = 25 * np.exp(-1.0 * x_val**2)
    plt.plot(x_val, gauss)
    plt.show()


if __name__ == "__main__"
main()
