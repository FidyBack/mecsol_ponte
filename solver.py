import numpy as np
from funcoesTermosol import *
from math import sqrt
from solve import solve
import itertools

# Nó
class Node:
    def __init__(self, number, x, y):
        self.x = x
        self.y = y
        self.liberty_degree = [(number*2)-1, number*2]

# Barra
class Bar:
    def __init__(self, node1, node2, modulus_of_elasticity, cross_section_area):
        self.liberty_degree = list(itertools.chain(*[node1.liberty_degree, node2.liberty_degree]))

        lenght = sqrt((node2.x - node1.x)**2 + (node2.y - node1.y)**2)
        s = (node2.y - node1.y) / lenght
        c = (node2.x - node1.x) / lenght
        MR = np.array ([[c**2, c*s, -c**2, -c*s],
                        [c*s, s**2, -c*s, -s**2],
                        [-c**2, -c*s, c**2, c*s],
                        [-c*s, -s**2, c*s, s**2]])
        
        self.rigidity_matrix = ((modulus_of_elasticity * cross_section_area) / lenght)*MR
        self.sin = s
        self.cos = c
        self.lenght = lenght
        self.array = np.array([-self.cos, - self.sin, self.cos, self.sin])
        self.modulus_of_elasticity = modulus_of_elasticity

    def calculate(self, displacements_vector):
        liberty_degree = [x - 1 for x in self.liberty_degree]
        displacements_vector = [displacements_vector[i] for i in liberty_degree]
        # return deformation, tension
        return (1/self.lenght) * self.array.dot(displacements_vector), (self.modulus_of_elasticity/self.lenght) * self.array.dot(displacements_vector)
    

class Solver:

    def __init__(self):
        self.solved = False

    def load(self, fille_name):
        data = importa(fille_name + ".xlsx")
        self.nodes_number = data[0]
        self.nodes_matrix = data[1]
        self.members_number = data[2]
        self.incidence_matrix = data[3]
        self.loads_number = data[4]
        self.loads_vector = data[5]
        self.restrictions_number = data[6]
        self.restrictions_vector = [int(x) for x in data[7]]
        self.solved = False

    def plot(self):
        plota(self.nodes_matrix, self.incidence_matrix)

    def solve(self):

        nodes = []
        for node in range(self.nodes_number):
            nodes.append(Node(node+1, self.nodes_matrix[0][node], self.nodes_matrix[1][node]))
        
        bars = []
        global_rigidity_matrix = np.zeros((self.nodes_number * 2, self.nodes_number * 2))
        
        
        for member in range(self.members_number):
            bar = Bar(nodes[(int(self.incidence_matrix[member][0]))-1], nodes[(int(self.incidence_matrix[member][1]))-1], self.incidence_matrix[member][2], self.incidence_matrix[member][3])
            bars.append(bar)
            for line in range(4):
                for column in range(4):
                    global_rigidity_matrix[bar.liberty_degree[line]-1][bar.liberty_degree[column]-1] += bar.rigidity_matrix[line][column]

        for l in range(len(global_rigidity_matrix)):
            print(global_rigidity_matrix[l])
        #Condições de Contorno
        contour_loads_vector = np.delete(self.loads_vector, self.restrictions_vector, 0)
        contour_global_rigidity_matrix = np.delete(global_rigidity_matrix, self.restrictions_vector, 0)
        contour_global_rigidity_matrix = np.delete(contour_global_rigidity_matrix, self.restrictions_vector, 1)
        
        contour_displacements_vector = solve(contour_global_rigidity_matrix, contour_loads_vector)

        displacements_vector = np.zeros(self.nodes_number * 2)
        index = np.arange(0, self.nodes_number*2, 1)
        index = np.delete(index, self.restrictions_vector, 0)

        for i in range(len(index)):
            displacements_vector[index[i]] = contour_displacements_vector[i]

        resultant_loads_vector = global_rigidity_matrix.dot(displacements_vector)

        tensions_vector = np.zeros(len(bars))
        deformations_vector = np.zeros(len(bars))

        for i in range(len(bars)):
            deformation, tension = bars[i].calculate(displacements_vector)
            tensions_vector[i] = tension
            deformations_vector[i] = deformation

        self.resultant_loads_vector = resultant_loads_vector
        self.deformations_vector = deformations_vector
        self.tensions_vector = tensions_vector
        self.displacements_vector= displacements_vector

        self.solved = True

    def write(self, fille_name):

        if self.solved:
            geraSaida(fille_name, self.resultant_loads_vector, self.displacements_vector, self.deformations_vector, [], self.tensions_vector)
        else:
            print("Solve first")

solver = Solver()
solver.load("entrada_quadrado")
solver.solve()
solver.write("saida_quadrado")