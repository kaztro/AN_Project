# -*- coding: UTF-8 -*-
 
'''
Cuadratura de Gauss-Legendre
Abraham Terán
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

