#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Este codigo aproxima la constante de Hubble con los segundos datos
experimentales. el codigo es muy similar a parte1.py para los datos iniciales
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
    datos = np.loadtxt(nombre, usecols=(80, 81, 82, 83))
    return datos


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


# main
nom = "data/DR9Q.dat"
datos = leer_archivo(nom)
i = datos[:, 0]
di = datos[:, 1]
z = datos[:, 2]
dz = datos[:, 3]
adivinanza = 500
# aprox via leastsq
aprox1 = aprox_leastsq(i, z, adivinanza)
print aprox1
# datos para graficar
# datos para graficar
i_min = np.amin(i)
i_max = np. amax(i)
i_aprox = np.linspace(i_min, i_max, 10)
z_aprox2 = y_aprox(i_aprox, aprox1[0])
# graficos
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(i, z, 'o')
ax1.plot(i_aprox, z_aprox2)
ax1.set_xlabel("d")
ax1.set_ylabel("v(d)")
plt.savefig("parte1.png")
plt.draw()
plt.show()