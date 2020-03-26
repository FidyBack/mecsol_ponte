# -*- coding: utf-8 -*-
import funcoesTermosol as ft
import numpy as np

tolerancia = 10**-64

def solve_line(rigidez, forcas, deslocamento, indice):
    a = np.concatenate((rigidez[indice][0:indice], rigidez[indice][indice+1:len(deslocamento)+1]))
    x = np.concatenate((deslocamento[0:indice], deslocamento[indice+1:len(deslocamento)+1]))
    b = forcas[indice][0]
    div = rigidez[indice][indice]
    sum = b
    for i in range(len(a)):
        sum -= a[i]*x[i]
    return sum / div

def solve(rigidez, forcas):
    n = len(forcas)
    deslocamento = np.zeros((n, 1))
      
    for i in range(n):
        deslocamento[i] =  solve_line(rigidez, forcas, deslocamento, i)

    dif_boll = [False]
    i = 0
    while(False in dif_boll):
        
        if(i > n-1):
            i = 0

        deslocamento_dif = deslocamento.copy()
        deslocamento[i] =  solve_line(rigidez, forcas, deslocamento, i)
        dif = deslocamento - deslocamento_dif
        dif_boll = np.logical_and(dif<  tolerancia, dif> -tolerancia)

        i += 1

    return deslocamento