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

# funciones estructurales


def leer_archivo(nombre):
    '''
    lee el archivo
    nombre debe ser un str
    '''
    datos = np.loadtxt(nombre)
    return datos
    pass


def aproximar_parametros(a):
    '''
    aproxima los parametos a. puede ser un escalar o un conjunto de escalares
    a optimizar
    '''
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


# main
nom = "data/hubble_original.dat"
datos = leer_archivo(nom)
d = datos[:, 0]
v = datos[:, 1]
adivinanza = 500
# aprox via leastsq
aprox1 = aprox_leastsq(d, v, adivinanza)
# aproximacion manual
aprox2 = aprox_manual(d, v, adivinanza)
print aprox2
# datos para graficar
d_aprox = np.linspace(d[0], d[-1], 10)
print aprox1[0]
v_aprox2 = y_aprox(d_aprox, aprox2)
# graficos
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(d, v, 'o')
ax1.plot(d_aprox, v_aprox2)
ax1.set_xlabel("d")
ax1.set_ylabel("v(d)")
plt.savefig("parte1.png")
plt.draw()
plt.show()
