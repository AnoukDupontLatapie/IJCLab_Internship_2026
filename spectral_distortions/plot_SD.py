import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## CONSTANTES##
h = 6.62607015e-34  # Planck constant
c = 299792458  # speed of light
k_B = 1.380649e-23  # Boltzmann constant
T = 2.725  # CMB temperature in Kelvin


def black_body(x):
    return 2*h/c**2 * nu**3/(np.exp(x)-1)


def temp_shift(x):
    return x*np.exp(x)/(np.exp(x)-1)**2


def y_dist(x):
    return temp_shift(x) * (x*(np.exp(x)+1)/(np.exp(x)-1) - 4)


def mu_dist(x):
    return temp_shift(x) * (0.4561 - 1/x)


nu = np.linspace(1, 1000, 1000) * 10**9  # GHz
x = h*nu/(k_B*T)

I_temp = 2*h*nu**3 * temp_shift(x) / c**2
I_y = 2*h*nu**3 * y_dist(x) / c**2
I_mu = 2*h*nu**3 * mu_dist(x) / c**2

I0 = 2*h/c**2 * (k_B*T/h)**3

plt.plot(x, 1/4 * I_temp / I0, linestyle='dashed', label='temperature shift')
plt.plot(x, 1/4 * I_y / I0, linestyle='dotted', label='y-distortion')
plt.plot(x, 1.401 * I_mu / I0, linestyle='dashdot', label='mu-distortion')
plt.plot(x, black_body(x)/I0, linestyle='solid', label='black body')
plt.xscale('log')
plt.xlabel(r'$x$', fontsize=14)
plt.ylabel(r'$\frac{I}{I0}$', fontsize=14)
plt.title('Spectral distortions scaled')
plt.legend()
plt.show()


# Redshift taken into account

def J(z):
    return np.exp(-(z/(2*10**6))**(5/2))


def J_mu(z):
    return 1 - np.exp(-((1+z)/(5.8*10**(-4)))**1.88)


def J_y(z):
    return (1 + ((1+z)/(6*10**4))**2.58)**(-1)


def green_func(nu, z):
    a = 2*h*nu**3 / c**2
    x = h*nu/(k_B*T)
    return a * (1.401 * mu_dist(x)/I0 * J_mu(z)*J(z) + 1/4 * y_dist(x)/I0 * J_y(z) + temp_shift(x)/(4*I0) * (1 - J(z)))


for z in np.linspace(2*10**6, 10**2, 100):
    plt.clf()
    plt.plot(nu, green_func(nu, z), label=f'z={z:.0e}')
    plt.plot(nu, green_func(nu, 2*10**6),
             linestyle='dotted', label=r'Black body')
    plt.plot(nu, green_func(nu, 3*10**5), linestyle='dotted',
             label=r'transition $\mu$ to y')
    plt.plot(nu, green_func(nu, 5*10**3),
             linestyle='dotted', label=r'y-distortion')
    plt.xscale('log')
    plt.xlabel(r'$\nu$ (GHz)', fontsize=14)
    plt.ylabel(r'$G_{th}$', fontsize=14)
    plt.legend()
    plt.pause(0.001)

plt.show()
