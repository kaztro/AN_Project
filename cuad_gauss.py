# -*- coding: UTF-8 -*-
 
'''
Cuadratura de Gauss-Legendre
Abraham Terán
'''
''' 
import math
from sympy import *
 
x = Symbol('x')
y = Symbol('y')
 
def rodrigues(n):   # Aqui 'n' es el grado del polinomio de Legendre
    y = (x**2 - 1)**n
    pol = diff(y,x,n)/(2**n * math.factorial(n)) # Fórmula de Rodrigues
    return pol
 
def main():
    func = raw_input("Ingrese funcion (y = ...) > ")
    A = float(input("Ingrese 'a' de rango  de integración (a<->b) > "))
    B = float(input("Ingrese 'b' de rango  de integración (a<->b) > "))
    n = input("Ingrese 'n' (n >= 2) > ")\
 
    file = open("data.py", "w")
    file.close()
    file = open("data.py", "a")
    file.write("import math\n")
    file.write("from sympy import *\n")
    file.write("x = Symbol('x')\n")
    file.write("y = Symbol('y')\n")
    file.write(func)
    file.close()\
 
    xIpots = solve(rodrigues(n), x) # Raices de polinomios de Legendre
    LePolD = diff(rodrigues(n))     # Derivada de polinomios de Legendre
    Cis = []
    DataFSum = [[],[]]\
 
    import data\
 
    for i in range(n):
        Cis.append(2/((1-xIpots[i]**2)*(LePolD.evalf(subs={x:xIpots[i]}))**2))
        DataFSum[0].append(data.y.evalf(subs={x:(B-A)*(xIpots[i]/2)+(A+B)/2}))
        DataFSum[1].append(Cis[i]*DataFSum[0][i])\
 
    cuad = ((B-A)/2)*sum(DataFSum[1])
    exac = (integrate(data.y,x).evalf(subs={x:B}))-(integrate(data.y,x).evalf(subs={x:A}))
    erro = abs((exac-cuad)*100)
    print "\nResultado de cuadratura: ...", cuad
    print "Resultado exacto : .........", exac
    print "Porcentaje de Error: .......", erro
    raw_input("\nENTER para Salir...")
 
print "\n\tCUADRATURA DE GAUSS\n"
main()
'''
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
'''
def switchN(f, a, b, C, x, n):
   return {
   	2: (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1])),
	3: (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2])),
	4: (b - a)/2 * (C[0] * f(x[0]) + C[1] * f(x[1]) + C[2] * f(x[2]) + C[3] * f(x[3])),
   }[n]
'''
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
