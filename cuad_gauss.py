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

#Grado 5
xin5 = np.array([-0.9061798459, -0.5384693101, 0,0.5384693101, 0.9061798459])
cin5 = np.array([0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705,0.2369268851])

#Grado 6
xin6 = np.array([-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142])
cin6 = np.array([0.1713244924, 0.3607615730, 0.4679139346, 0.4679139346,  0.3607615730, 0.1713244924])

#Grado 7
xin7 = np.array([-0.9491079123, -0.7415311856, -0.4058451514, 0, 0.4058451514, 0.7415311856, 0.9491079123])
cin7 = np.array([0.1294849662, 0.2797053915, 0.3818300505, 0.4179591837, 0.3818300505, 0.2797053915, 0.1294849662])

#Grado 8
xin8 = np.array([-0.9602898565, -0.7966664774, -0.5255324099, -0.1834346425, 0.1834346425, 0.5255324099, 0.7966664774, 0.9602898565])
cin8 = np.array([0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834, 0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363])

#Grado 9
xin9 = np.array([-0.9681602395, -0.8360311073, -0.6133714327, -0.3242534234, 0, 0.3242534234, 0.6133714327, 0.8360311073, 0.9681602395])
cin9 = np.array([0.0812743883, 0.1806481607, 0.2606106964, 0.3123470770, 0.3302393550, 0.3123470770, 0.2606106964, 0.1806481607,0.0812743883])

#Grado 10
xin10 = np.array([-0.9739065285, -0.8650633667, -0.6794095683, -0.4333953941, -0.1488743390, 0.1488743390, 0.4333953941, 0.6794095683, 0.8650633667, 0.9739065285])
cin10 = np.array([0.0666713443, 0.1494513492, 0.2190863625, 0.2692667193, 0.2955242247, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513492, 0.0666713443])

#gauss_legendre(function, a, b, Root, Coefficient, n)
def gauss_legendre(f, a, b, R, C, n):
    x = np.zeros(n)
    for i in range(n):
        x[i] = (b + a)/2 + (b - a)/2 * R[i]

    if (n == 10): result = (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2]) + C[3] * f(x[3]) + C[4] * f(x[4]) + C[5] * f(x[5]) + C[6] * f(x[6]) + C[7] * f(x[7]) + C[8] * f(x[8])  + C[9] * f(x[9]))
    if (n == 9): result = (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2]) + C[3] * f(x[3]) + C[4] * f(x[4]) + C[5] * f(x[5]) + C[6] * f(x[6]) + C[7] * f(x[7])  + C[8] * f(x[8]))
    if (n == 8): result = (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2]) + C[3] * f(x[3]) + C[4] * f(x[4]) + C[5] * f(x[5]) + C[6] * f(x[6]) + C[7] * f(x[7]))
    if (n == 7): result = (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2]) + C[3] * f(x[3]) + C[4] * f(x[4]) + C[5] * f(x[5]) + C[6] * f(x[6]))
    if (n == 6): result = (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2]) + C[3] * f(x[3]) + C[4] * f(x[4]) + C[5] * f(x[5]))
    if (n == 5): result = (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2]) + C[3] * f(x[3]) + C[4] * f(x[4]))
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
Gau5 = gauss_legendre(f, a, b, xin5, cin5, 5)
Gau6 = gauss_legendre(f, a, b, xin6, cin6, 6)
Gau7 = gauss_legendre(f, a, b, xin7, cin7, 7)
Gau8 = gauss_legendre(f, a, b, xin8, cin8, 8)
Gau9 = gauss_legendre(f, a, b, xin9, cin9, 9)
Gau10 = gauss_legendre(f, a, b, xin10, cin10, 10)

print("Gaussian integral n=2: ", Gau2)
print("Gaussian integral n=3: ", Gau3)
print("Gaussian integral n=4: ", Gau4)
print("Gaussian integral n=5: ", Gau5)
print("Gaussian integral n=6: ", Gau6)
print("Gaussian integral n=7: ", Gau7)
print("Gaussian integral n=8: ", Gau8)
print("Gaussian integral n=9: ", Gau9)
print("Gaussian integral n=10: ", Gau10)

#f = lambda x: x**4
#f = lambda x: (4*x + 4) / (np.sqrt( (2*x + 2)**2 + 9 ))
