# -*- coding: utf-8 -*-
"""
Assignment 5 problem 1
Computes the integral I=\int_0^1 \sin^2(\sqrt{100x}) dx
using the adaptive trapezoidal rule method. Part (a) uses
the adaptive trapezoidal rule (in commments), and for (b) the modified
program uses Romberg integration. The output in both 
cases is the number of slices, estimate of the integral, and the
estimate of the error on the integral.
"""
from math import sin,sqrt,fabs
from numpy import array
def f(x):
    return sin(sqrt(100.0*x))*sin(sqrt(100.0*x))
a=0.0
b=1.0
# N=1
epsilon=1.0e-6
# h=(b-a)/N
#I=h*(0.5*f(a)+0.5*f(b))       # initial estimate of the integral for one slice
err=1.0                     # to be sure the while loop evaluates at least once
# I_i=I
#print("For 1 slice the integral is approximately {0:.5f}".format(I))
#while fabs(err)>epsilon:
#    N*=2
#    h=(b-a)/N                          # PART (a) ADAPTIVE TRAPEZOIDAL METHOD 
#    s=0.0
#    for k in range(1,N,2):
#        s+=f(a+k*h)
#    I_i=0.5*I+h*s         # I_i estimate of the integral
#    err=(1.0/3.0)*(I_i-I)       # new error
#    print("For {0} slices the integral is approximately {1:.5f} with approximate error {2:.5e}".format(N,I,err))
#    I=I_i

#############################################################################
# PART (b) ROMBERG INTEGRATION
n=6                 # increase this to increase the number of levels of recursion
r = array([[0]*(n+1)]*(n+1),float)
h = b-a
r[0,0] =0.5*h*(f(a)+f(b))
N = 1
for i in range(1,n+1):
    h = 0.5*h
    sum = 0.0
    N = 2*N
    for k in range(1, N, 2):
        sum = sum +f(a+k*h)        
    r[i,0] = 0.5*r[i-1,0]+sum*h
    M = 1                               
    for j in range(1,i+1):
        M = 4*M                     # this is 4^m in the book
        r[i,j] = r[i,j-1] + (r[i,j-1]-r[i-1,j-1])/(M-1)
    err=1/(M-1)*(r[i,j-1]-r[i-1,j-1])
    if fabs(err)<epsilon: 
       # print("Required accuracy reached at {0} slices".format(N))
        break
print(r)            # prints out the triangular array of values (ala pg 161)
print("Romberg integration converges to the required accuracy at {0} slices".format(N))

#############################################################################
# OUTPUT OF PROGRAM FROM PART (a)
# For 1 slice the integral is approximately 0.14798
# For 2 slices the integral is approximately 0.14798 with approximate error 5.90841e-02
# For 4 slices the integral is approximately 0.32523 with approximate error 6.23503e-02
# For 8 slices the integral is approximately 0.51228 with approximate error -3.64285e-02
# For 16 slices the integral is approximately 0.40300 with approximate error 9.03531e-03
# For 32 slices the integral is approximately 0.43010 with approximate error 6.10377e-03
# For 64 slices the integral is approximately 0.44841 with approximate error 1.83276e-03
# For 128 slices the integral is approximately 0.45391 with approximate error 4.78524e-04
# For 256 slices the integral is approximately 0.45535 with approximate error 1.20921e-04
# For 512 slices the integral is approximately 0.45571 with approximate error 3.03111e-05
# For 1024 slices the integral is approximately 0.45580 with approximate error 7.58283e-06
# For 2048 slices the integral is approximately 0.45582 with approximate error 1.89602e-06
# For 4096 slices the integral is approximately 0.45583 with approximate error 4.74026e-07
#
# Adaptive trapezoidal integration converges to the required accuracy at 4096 slices
##############################################################################
# OUTOUT OF PROGRAM FROM PART (b)
#[[ 0.14797948  0.          0.          0.          0.          0.          0.        ]
# [ 0.32523191  0.38431605  0.          0.          0.          0.          0.        ]
# [ 0.51228285  0.57463317  0.58732097  0.          0.          0.          0.        ]
# [ 0.40299745  0.36656898  0.35269804  0.34897386  0.          0.          0.        ]
# [ 0.43010337  0.43913868  0.44397666  0.44542552  0.44580376  0.          0.        ]
# [ 0.44841467  0.45451843  0.45554375  0.45572735  0.45576775  0.45577749  0.        ]
# [ 0.45391293  0.45574569  0.4558275   0.45583201  0.45583242  0.45583248  0.45583249]]
# Romberg integration converges to the required accuracy at 64 slices
