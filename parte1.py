#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Este codigo aproxima la constante de Hubble con los primeros datos
experimentales. Se espera que no sea suficientemente cercano al valor
real
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy import optimize as opt

np.random.seed(8624617)
# funciones estructurales


def leer_archivo(nombre):
    '''
    lee el archivo
    nombre debe ser un str
    '''
    datos = np.loadtxt(nombre)
    return datos
    pass


def y_aprox(x, param):
    '''
    funcion lineal para minimizar
    se define para tener un caso mas general en el que se tenga una funcion
    mucho mas compleja por minimizar
    '''
    a = param
    b = 0
    y = a * x + b
    return y


def residuos(param, x, y):
    c = param
    y_ap = y_aprox(x, param)
    err = y - y_ap
    return err


def aprox_leastsq(d, v, adivinanza):
    aprox = leastsq(residuos, adivinanza, args=(d, v))
    return aprox


def aprox_manual(x, y, beta_0):
    '''
    se buscan los ceros de la derivada de chi cuadrado. es equivalente a usar
    la funcion leastsq.
    x e y son los arreglos de datos experimentales. beta es el parametro a
    minimizar. este caso es particular para una relacion lineal entre x e y de
    la forma y = m * x + n, con n = 0. m corresponde a beta
    '''
    derivada_chi_cuadrado = lambda beta: np.sum(x * y - x ** 2 * beta)
    beta = opt.newton(derivada_chi_cuadrado, beta_0)
    return beta


def intervalo_confianza(muestra_x, muestra_y, porcentaje):
    '''
    busca intervalo de confianza. genera una muestra aleatoria en base a los
    datos experimentales. a partir de cada muestra obtiene el valor de H_0.
    corresponde al metodo de bootstrapt. porcentaje refiere al porcentaje del
    intervalo de confianza que se busca.
    '''
    if (porcentaje < 0) or (porcentaje > 100):
        print str('porcentaje no valido')
    N = len(muestra_x)
    N_iteracion = int(N ** 2)
    promedios = np.zeros(N_iteracion)
    muestra_x, muestra_y
    for i in range(N_iteracion):
        azar = np.random.randint(low=0, high=N, size=N)
        x_i = muestra_x[azar]
        y_i = muestra_y[azar]
        aprox = biseccion(x_i, y_i)
        promedios[i] = aprox
    promedios = np.sort(promedios)
    minim = ((100 - porcentaje) /2) * 0.01
    maxim = 1 - (minim)
    lim_min = promedios[int(N_iteracion * minim)]
    lim_max = promedios[int(N_iteracion * maxim)]
    return lim_min, lim_max
    pass


def biseccion(d, v, adivinanza=500):
    aprox1 = aprox_manual(d, v, adivinanza)
    aprox2 = aprox_manual(v, d, 1. / adivinanza)
    ap1 = aprox1
    ap2 = 1 / aprox2
    a = (ap1 * ap2 - 1 + np.sqrt((1 + ap1 ** 2) * (1 + ap2 ** 2))) / (ap1 + ap2)
    return a
    pass


# main
nom = "data/hubble_original.dat"
datos = leer_archivo(nom)
d = datos[:, 0]
v = datos[:, 1]
adivinanza = 500
# aprox via leastsq
aprox1 = aprox_leastsq(d, v, adivinanza)
# aprox manual + biseccion
H_0 = biseccion(d, v, adivinanza)
print H_0
# datos para graficar
d_aprox = np.linspace(d[0], d[-1], 10)
print aprox1[0]
v_aprox2 = y_aprox(d_aprox, H_0)
# intervalo de confianza
intervalo = intervalo_confianza(d, v, 95)
print intervalo[0]
print intervalo[1]
# graficos
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(d, v, 'o')
ax1.plot(d_aprox, v_aprox2)
ax1.set_xlabel("Distancia")
ax1.set_ylabel("Velocidad")
plt.savefig("parte1.png")
plt.draw()
plt.show()
