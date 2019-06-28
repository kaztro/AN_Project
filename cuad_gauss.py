# -*- coding: UTF-8 -*-
 
import numpy as np

#Grado 2
xin2 = np.array([-0.57735026918963, 0.57735026918963])
cin2 = np.array([1.00000000000000, 1.00000000000000])

#Grado 3
xin3 = np.array([-0.77459666924148, 0.00000000000000, 0.77459666924148])
cin3 = np.array([0.55555555555555, 0.88888888888888, 0.55555555555555])

#Grado 4
xin4 = np.array([-0.86113631159405, -0.33998104358486, 0.33998104358486, 0.86113631159405])
cin4 = np.array([0.34785484513745, 0.65214515486255, 0.65214515486255, 0.34785484513745])

#gauss_legendre(function, a, b, Root, Coefficient, n)
def gauss_legendre(f, a, b, R, C, n):
    x = np.zeros(n)
    for i in range(n):
        x[i] = (b + a)/2 + (b - a)/2 * R[i]
     
    if (n == 4): result = (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2]) + C[3] * f(x[3]))
    if (n == 3): result = (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2]))
    if (n == 2): result = (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]))	    
	
    return result

#F function
f = lambda x: np.exp(-(x**2)/2)
#Interval
a = -1.0; b = 1.0

Gau2 = gauss_legendre(f, a, b, xin2, cin2, 2)
Gau3 = gauss_legendre(f, a, b, xin3, cin3, 3)
Gau4 = gauss_legendre(f, a, b, xin4, cin4, 4)

print("Gaussian integral n=2: ", Gau2)
print("Gaussian integral n=3: ", Gau3)
print("Gaussian integral n=4: ", Gau4)
