import numpy as np
from funcoesTermosol import *
from math import sqrt
from solve import solve
import itertools


# Funções de entrada
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xlsx')

# Nó
class Node:
    def __init__(self, number, x, y):
        self.number = number
        self.x = x
        self.y = y
        self.lib = [(number*2)-1, number*2]

# Barra
class Bar:
    def __init__(self, number, no1, no2, E, A):
        self.number = number
        self.no1 = no1
        self.no2 = no2
        self.E = E
        self.A = A
        self.L = sqrt((no2.x - no1.x)**2 + (no2.y - no1.y)**2)
        s = (no2.y - no1.y)/self.L
        c = (no2.x - no1.x)/self.L
        MR = np.array ([[c**2, c*s, -c**2, -c*s],
                        [c*s, s**2, -c*s, -s**2],
                        [-c**2, -c*s, c**2, c*s],
                        [-c*s, -s**2, c*s, s**2]])
        self.mat_rig = ((self.E*self.A)/self.L)*MR

# Calcula os vetores de deslocamento nodal de acordo com a quatidade de elementos passados
def vet_des_nod(nn, N, nm, Inc, nc, F, nr, R):
    # Variáveis
    nodes = []
    bars = []

    for i in range(nn):
        nodes.append(Node(i+1, N[0][i], N[1][i]))
    for i in range(nm):
        bars.append(Bar(i+1, nodes[(int(Inc[i][0]))-1], nodes[(int(Inc[i][1]))-1], Inc[i][2], Inc[i][3]))

    # Matriz Rigidez Global
    mrg = np.zeros((nm*2, nm*2))

    for i in range(len(bars)):
        mrg[bars[i].no1.lib[0]-1:bars[i].no1.lib[1] , bars[i].no1.lib[0]-1:bars[i].no1.lib[1]] += bars[i].mat_rig[0:2, 0:2]
        mrg[bars[i].no1.lib[0]-1:bars[i].no1.lib[1] , bars[i].no2.lib[0]-1:bars[i].no2.lib[1]] += bars[i].mat_rig[0:2, 2:4]
        mrg[bars[i].no2.lib[0]-1:bars[i].no2.lib[1] , bars[i].no1.lib[0]-1:bars[i].no1.lib[1]] += bars[i].mat_rig[2:4, 0:2]
        mrg[bars[i].no2.lib[0]-1:bars[i].no2.lib[1] , bars[i].no2.lib[0]-1:bars[i].no2.lib[1]] += bars[i].mat_rig[2:4, 2:4]

    # Condições de Contorno
    Fcc = np.delete(F, R, 0)
    MRlin = np.delete(mrg, R, 0)
    MRcc = np.delete(MRlin, R, 1)

    # Inversão e Multiplicação
    U = solve(MRcc, Fcc)
    return(U, mrg)

desloc = vet_des_nod(nn, N, nm, Inc, nc, F, nr, R)

print(desloc[0])