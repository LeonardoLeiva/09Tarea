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
    return datos * 3.631


def y_aprox(x, param):
    '''
    funcion lineal para minimizar
    se define para tener un caso mas general en el que se tenga una funcion
    mucho mas compleja por minimizar
    a es pendiente, b es coef de posicion
    '''
    a, b = param
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


def intervalo_confianza(muestra_x, muestra_y, err_x, err_y, porcentaje):
    '''
    busca intervalo de confianza. genera una muestra aleatoria en base a los
    datos experimentales. a partir de cada muestra obtiene el valor de la
    constante apropiada para modelo lineal.
    corresponde al metodo de monte carlo. porcentaje refiere al porcentaje del
    intervalo de confianza que se busca.
    '''
    N = len(muestra_x)
    Nmc = 10000
    promedios_a = np.zeros(Nmc)
    promedios_b = np.zeros(Nmc)
    for i in range(Nmc):
        r = np.random.normal(0, 1, size=len(muestra_x))
        x_i = muestra_x + err_x * r
        y_i = muestra_y + err_y * r
        a_i, b_i = biseccion(x_i, y_i)
        promedios_a[i-1] = a_i
        promedios_b[i-1] = b_i
    promedios_a = np.sort(promedios_a)
    promedios_b = np.sort(promedios_b)
    minim = ((100 - porcentaje) /2) * 0.01
    maxim = 1 - (minim)
    lim_min_a = promedios_a[int(Nmc * minim)]
    lim_max_a = promedios_a[int(Nmc * maxim)]
    lim_min_b = promedios_b[int(Nmc * minim)]
    lim_max_b = promedios_b[int(Nmc * maxim)]
    return lim_min_a, lim_max_a, lim_min_b, lim_max_b
    pass


def biseccion(x, y):
    a1, b1 = np.polyfit(x, y, 1)
    a_0, b_0 = np.polyfit(y, x, 1)
    a2 = 1 / a_0
    b2 = - b_0 / a_0
    x_0 = (b2 - b1) / (a1 - a2)
    y_0 = (a1 * b2 - b1 * a2) / (a1 - a2)
    a = (a1 * a2 - 1 + np.sqrt((1 + a1 ** 2) * (1 + a2 ** 2))) / (a1 + a2)
    b = (y_0 - a * x_0)
    return a, b


# main
nom = "data/DR9Q.dat"
datos = leer_archivo(nom)
i = datos[:, 0]
di = datos[:, 1]
z = datos[:, 2]
dz = datos[:, 3]
adivinanza = 1
# aprox via leastsq
aprox1 = biseccion(i, z)
print aprox1
# datos para graficar
i_min = np.amin(i)
i_max = np. amax(i)
i_aprox = np.linspace(i_min, i_max, 10)
z_aprox2 = y_aprox(i_aprox, aprox1)
# intervalo de confianza
intervalo = intervalo_confianza(i, z, di, dz, 95)
print intervalo
# graficos
fig = plt.figure()
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(i, z, 'o')
ax1.plot(i_aprox, z_aprox2)
ax1.set_xlabel("Flujo de banda i")
ax1.set_ylabel("Flujo de banda z")
plt.savefig("parte1.png")
plt.draw()
plt.show()
