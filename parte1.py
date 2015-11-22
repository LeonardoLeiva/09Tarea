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


# main
nom = "data/hubble_original.dat"
datos = leer_archivo(nom)
d = datos[:, 0]
v = datos[:, 1]
adivinanza = 500
aprox = leastsq(residuos, adivinanza, args=(d, v))
# datos para graficar
d_aprox = np.linspace(d[0], d[-1], 10)
H_0 = aprox[0]
print H_0
v_aprox = y_aprox(d_aprox, H_0)
# graficos
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(d, v, 'o')
ax1.plot(d_aprox, v_aprox)
ax1.set_xlabel("d")
ax1.set_ylabel("v(d)")
plt.savefig("parte1.png")
plt.draw()
plt.show()
