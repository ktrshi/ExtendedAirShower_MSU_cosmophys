import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from iminuit import Minuit
from scipy import integrate
import mplcatppuccin
from mplcatppuccin.palette import load_color
from scipy.optimize import curve_fit

mpl.style.use("latte")


def func(x, k):
    return c[k] * X ** a[k]


x = [1, 10]
Fe = [5.836e7, 6.494e8]
Fe_er = [6.87e6, 1.043e8]
P = [8.504e7, 6.965e8]
P_er = [1.439e7, 2.056e8]

a = np.zeros(2)
b = np.zeros(2)
c = np.zeros(2)
var_list = [Fe, P]
for k, var in enumerate(var_list):
    log_x = np.log(x)
    log_y = np.log(var)
    a[k], b[k] = np.polyfit(log_y, log_x, 1)
    c[k] = math.exp(b[k])
X = np.arange(5e7, 7e8, 5000)
Y = c[0] * X ** a[0]
plt.figure()
# plt.scatter(X, Y)
# color = load_color("mocha", "peach")
plt.errorbar(Fe, x, xerr=Fe_er, fmt='o', capsize=3,
             label='Fe\nC = {:.3}\n'.format(c[0]) + r'$\beta$' + ' = {:.3}'.format(a[0]))
plt.plot(X, func(X, 0), label='Fe')
plt.errorbar(P, x, xerr=P_er, fmt='o', capsize=3,
             label='P\nC = {:.3}\n'.format(c[1]) + r'$\beta$' + ' = {:.3}'.format(a[1]))
plt.plot(X, func(X, 1), label='P')
plt.ylabel('E, ПэВ')
plt.xlabel(r'$\int$ 300')
plt.title(r'$E_{0}(Q_{ch}300)$')
plt.legend()
plt.savefig('final1.jpg', dpi=200, bbox_inches='tight')

# X = np.arange(1, 1e6, 2e2)
# Y =
# plt.figure()
# plt.plot(X, Y)
# plt.show()
