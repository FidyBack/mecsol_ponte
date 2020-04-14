import numpy as np
from funcoesTermosol import *
from math import sqrt
from solve import solve
import itertools

# Nó
class Node:
    def __init__(self, number, x, y):
        self.number = number
        self.x = x
        self.y = y
        self.liberty_degree = [(number*2)-1, number*2]

# Barra
class Member:
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
        self.lenght = lenght
        self.array = np.array([-c, -s, c, s])
        self.modulus_of_elasticity = modulus_of_elasticity
        self.cross_section_area = cross_section_area

    def calculate(self, displacements_vector):
        liberty_degree = [x - 1 for x in self.liberty_degree]
        displacements_vector = [displacements_vector[i] for i in liberty_degree]
        # return deformation, tension, interna forces
        deformation = (1/self.lenght) * self.array.dot(displacements_vector)
        tension = self.modulus_of_elasticity * deformation
        internal_forces = tension * self.cross_section_area
        return deformation, tension, internal_forces
    

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

        self.nodes = []
        for node in range(self.nodes_number):
            self.nodes.append(Node(node+1, self.nodes_matrix[0][node], self.nodes_matrix[1][node]))
        
        self.members = []
        for member in range(self.members_number):
            member = Member(self.nodes[(int(self.incidence_matrix[member][0]))-1], self.nodes[(int(self.incidence_matrix[member][1]))-1], self.incidence_matrix[member][2], self.incidence_matrix[member][3])
            self.members.append(member)

    def plot(self):
        plota(self.nodes_matrix, self.incidence_matrix)

    def solve(self):

        global_rigidity_matrix = np.zeros((self.nodes_number * 2, self.nodes_number * 2))
        for member in self.members:
            for line in range(4):
                for column in range(4):
                    global_rigidity_matrix[member.liberty_degree[line]-1][member.liberty_degree[column]-1] += member.rigidity_matrix[line][column]

        #Countour Conditions
        contour_loads_vector = np.delete(self.loads_vector, self.restrictions_vector, 0)
        contour_global_rigidity_matrix = np.delete(global_rigidity_matrix, self.restrictions_vector, 0)
        contour_global_rigidity_matrix = np.delete(contour_global_rigidity_matrix, self.restrictions_vector, 1)
        
        #Solve for U. Numeric solve method developed, and on solve.py file. Numpy's is faster
        contour_displacements_vector = np.linalg.solve(contour_global_rigidity_matrix, contour_loads_vector)

        #Global U
        displacements_vector = np.zeros((self.nodes_number * 2, 1))
        index = np.arange(0, self.nodes_number*2, 1)
        index = np.delete(index, self.restrictions_vector, 0)
        for i in range(len(index)):
            displacements_vector[index[i]] = contour_displacements_vector[i]


        for i in range(self.nodes_number): # Deslocamentos multiplicados por 100 para visualização
            self.nodes_matrix[0][i] += displacements_vector[i*2]
            self.nodes_matrix[1][i] += displacements_vector[1 + i*2]
        #Loads vector = Kg * Ug
        resultant_loads_vector = global_rigidity_matrix.dot(displacements_vector)
        target_loads_vector = np.delete(resultant_loads_vector, index, 0)
        

        tensions_vector = np.zeros((self.members_number, 1))
        deformations_vector = np.zeros((self.members_number, 1))
        internal_forces_vector = np.zeros((self.members_number, 1))

        #Tension and Deformation 
        for i in range(self.members_number):
            deformation, tension, internal_forces = self.members[i].calculate(displacements_vector)
            tensions_vector[i] = tension
            deformations_vector[i] = deformation
            internal_forces_vector[i] = internal_forces

        self.resultant_loads_vector = resultant_loads_vector
        self.target_loads_vector = target_loads_vector
        self.deformations_vector = deformations_vector
        self.tensions_vector = tensions_vector
        self.internal_forces_vector = internal_forces_vector
        self.displacements_vector= displacements_vector

        self.solved = True

    def write(self, fille_name):

        if self.solved:
            geraSaida(fille_name, self.target_loads_vector, self.displacements_vector, self.deformations_vector, self.internal_forces_vector, self.tensions_vector)
        else:
            print("Solve first")

solver = Solver()

solver.load("entradas/entrada_ponte")
solver.plot()
solver.solve()
solver.plot()
solver.write("saidas/saida_ponte")