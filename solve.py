# -*- coding: utf-8 -*-
# sudo nc -l -p 888 -vvv
import funcoesTermosol as ft
import numpy as np

tolerancia = 10**-32

def solve_line(rigidez, forcas, deslocamento, indice):
    a = np.concatenate((rigidez[indice][0:indice], rigidez[indice][indice+1:len(deslocamento)+1]))
    x = np.concatenate((deslocamento[0:indice], deslocamento[indice+1:len(deslocamento)+1]))
    b = forcas[indice]
    div = rigidez[indice][indice]
    sum = b
    for i in range(len(a)):
        sum -= a[i]*x[i]
    return sum / div

def solve(rigidez, forcas):
    deslocamento = np.zeros(len(forcas))
    dif_boll = [False]   

    while(False in dif_boll):
        deslocamento_dif = deslocamento.copy()

        for i in range(len(forcas)):
            deslocamento[i] =  solve_line(rigidez, forcas, deslocamento, i)


        dif = deslocamento - deslocamento_dif
        dif_boll = np.logical_and(dif<  tolerancia, dif> -tolerancia)


    return deslocamento

rigidez = np.array([[1.59, -0.40,-0.54], 
                    [-0.40, 1.70, 0.40], 
                    [-0.54, 0.40, 0.54]])
rigidez *= 10**8
forcas = np.array([0, 150, -100])

print(solve(rigidez, forcas))