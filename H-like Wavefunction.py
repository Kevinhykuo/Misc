import sys
import numpy as np
import matplotlib.pyplot as plt
import math
from math import sqrt
from numpy import sin, cos, pi, exp
from scipy.special import eval_laguerre, eval_legendre, factorial
from scipy.misc import derivative
from mpl_toolkits.mplot3d import Axes3D

# The quantum numbers and the atomic mass
n = 2    #Principle quantum number
l = 0    #Orbital angular momentum quantum number
m = 0    #Magnetic quantum number
z = 1    #Atomic number

# Some constants
epsilon0 = 8.854187817*10**(-12)    #Permittivity of free space
hbar = 1.054571726*10**(-34)        #Plank's constant
e = 1.602176487*10**(-19)           #Electron charge
mp = 1.672621637*10**(-27)          #Mass of proton
me = 9.10938215*10**(-31)           #Mass of electron

# Some constant calculations
rm = mp*me / (mp+me)                        #Reduced mass
a0 = 4*pi*epsilon0*hbar**2 / (z*rm*e**2)   #Bohr radius
# a0 = 1

c1 = 1/sqrt(2*pi)
c2 = sqrt( (2*l+1)*factorial(l-abs(m)) / 2*factorial(l+abs(m)) )
c3 = sqrt( 4*factorial(n-l-1)*z**3 / (factorial(n+l)**3 * n**4 * a0**3) )
C = c1*c2*c3

a = np.linspace(0, 2*pi, 100)
p = np.linspace(0, pi, 100)
r = np.linspace(0, a0*10, 100)

# ===============================
# Azimuthal angle component of WF
# ===============================
def A(a):
    return( exp(1j*m*a) )

# Probability amplitude
def Apa(a):
    apa = abs(A(a))**2
    return(apa)

# ===========================
# Polar angle component of WF
# ===========================
def P1(p):
    p1 = sin(p)**abs(m)
    return(p1)

def P2(p):
    p2 = eval_legendre(l, cos(p))
    return(p2)

def DP2(p):
    dp2 = derivative(P2(p), p, abs(m))
    return(dp2)

def P(p):
    p = P1(p) * DP2(p)
    return(p)

# Probability amplitude
def Ppa(p):
    ppa = p**2
    return(ppa)

# ======================
# Radial component of WF
# ======================
def R1(r):
    return( ((2*z*r)/(n*a0))**l )

def R2(r):
    return( exp(-z*r/(n*a0)) )

def R3(r):
    #return( 0.5*((2*z*r/(n*a0))**2 - 4*(2*z*r/(n*a0)) + 2) )
    return( eval_laguerre(n+l, 2*z*r/(n*a0)) )

def DR3(r):
    dr3 = derivative(R3, 2*z*r/(n*a0), 2*l+1)
    return(dr3)

def R(r):
    R = R1(r)*R2(r)*DR3(r)
    return(R)

# # Probability amplitude
def Rpa(r):
    rpa = R(r)**2
    return(rpa)

# =========================================
# The wavefunction for a hydrogen-like atom
# =========================================
def WF(r, p, a):
    wf = R(r)*A(a)*P(p)*C
    return(wf)

# Probability amplitude
def PA(r, p, a):
    pa = Rpa(r)*Apa(a)*Ppa(p)*(C**2)
    return(pa)


#=================
# Plotting time!!!
#=================
y = a0*(r*R(r))**2
x = r/a0
# plt.plot(r, Rpa(r))
plt.plot(x, y)

# ax = plt.subplot(projection='polar')
# fig = plt.figure()
# ax2 = plt.add_subplot(111, projection='3d')
# ax.plot(p, Ppa(p))
# ax.plot(a, Apa(a))
#
plt.show()
