# Assignment_5
PHYS 50733 - Assignment 5 Solutions

problem1.py computes the integral I=\int_0^1 \sin^2(\sqrt{100x}) dx using the adaptive
trapezoidal method (in comments), and using Romberg integration. The output in both cases is the number of slices,
estimate of the integral, and the estimate of the error on the integral.

problem2.py computes the value of f(x)=1+tanh(2x),
then uses the central difference to calculate 
the derivative of the function in the range 
-2<=x<-2. The analytic formula for the derivative
of tanh(x) is sech^2(x), where sech(x) is the hyperbolic secant.
Thus, f'(x)=2sech(2x). We plot this formula for the aforementioned 
range of x-values on the same plot as the numerical answer. 

problem3.py computes the total energy per unit area
radiated by a blackbody by computing the integral
I=\int_0^\infty \frac{x^3}{\exp{x}-1} dx multiplied 
by a constant. Analytically this integral may in fact be evaluated exactly, 
using the Riemann zeta and Gamma functions. Using the value for the integral
we can estimate the Stefan Boltzmann constant, which has the approximate 
value \sigma ~ 5.670373(21)e-8 W m^-2 K^-4. We use Gaussian quadrature to 
compute the integral.

problem3_analytic.pdf derives the integral, computes it's exact value, and computes Stefan-Boltzmann constant.
