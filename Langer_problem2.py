# -*- coding: utf-8 -*-
"""
Assignment 5 problem 2
Computes the value of f(x)=1+tanh(2x),
then uses the central difference to calculate 
the derivative of the function in the range 
-2<=x<-2. The analytic formula for the derivative
of tanh(x) is sech^2(x), where sech(x) is the hyperbolic secant.
Thus, f'(x)=2sech(2x). We plot this formula for the aforementioned 
range of x-values on the same plot as the numerical answer. 
"""
from math import tanh,cosh
from pylab import plot,xlabel,ylabel,show,legend
def f(x):
    return 1.0+tanh(2.0*x)
def deriv_f(x):                 # for analytic solution
    return 2.0/(cosh(2*x)*cosh(2*x))
def Df(x):
    h=1.0e-5                # using h value recommended by book
    return (f(x+h/2.0)-f(x-h/2.0))/h
x=-2.0
X=[]
Y1=[]
Y2=[]
while x<2.0:
    X.append(x)
    Y1.append(deriv_f(x))           # analytic
    Y2.append(Df(x))                # numerical (central difference)
    x+=0.05
plot(X,Y1,"k--",label="analytic method")
plot(X,Y2,"go", label="central difference method")
xlabel("x")
ylabel("df/dx")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
show()