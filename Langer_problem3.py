# -*- coding: utf-8 -*-
"""
Assignment 5 problem 3
Computes the total energy per unit area
radiated by a blackbody by computing the integral
I=\int_0^\infty \frac{x^3}{\exp{x}-1} dx multiplied 
by a constant. Analytically this integral may in fact be evaluated exactly, 
using the Riemann zeta and Gamma functions. Using the value for the integral
we can estimate the Stefan Boltzmann constant, which has the approximate 
value \sigma ~ 5.670373(21)e-8 W m^-2 K^-4. We use Gaussian quadrature to 
compute the integral.
"""
from math import exp,fabs
from numpy import ones,copy,cos,tan,pi,linspace
def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

def f(z): 
    if z<1.0e-4:
        return z**2+7*z**3/2.0+91*z**4/12+13*z**5+13859*z**6/720+2309*z**7/90
    else:
        return z**3.0/((1.0-z)**5.0*(exp(z/(1-z))-1))
N=30
a=0.0
b=1.0
x,w=gaussxwab(N,a,b)
s=0.0
for k in range(N):
    s+=w[k]*f(x[k])
print("Value of integral is {0:.3f}".format(s))
#########################################################################
# OUTPUT OF PROGRAM (value of integral)
# Value of integral is 6.493947736403491
# 
