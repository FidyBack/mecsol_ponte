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
        self.liberty_degree = [(number*2), 1 + number*2]

# memberra
class Member:
    def __init__(self, node1, node2, modulus_of_elasticity, cross_section_area, density):
        self.liberty_degree = list(itertools.chain(*[node1.liberty_degree, node2.liberty_degree]))

        lenght = sqrt((node2.x - node1.x)**2 + (node2.y - node1.y)**2)
        s = (node2.y - node1.y) / lenght
        c = (node2.x - node1.x) / lenght
        MR = np.array ([[c**2, c*s, -c**2, -c*s],
                        [c*s, s**2, -c*s, -s**2],
                        [-c**2, -c*s, c**2, c*s],
                        [-c*s, -s**2, c*s, s**2]])
        self.node_1 = node1
        self.node_2 = node2
        self.sin = s
        self.cos = c
        self.rigidity_matrix = ((modulus_of_elasticity * cross_section_area) / lenght)*MR
        self.lenght = lenght
        self.array = np.array([-c, -s, c, s])
        self.modulus_of_elasticity = modulus_of_elasticity
        self.cross_section_area = cross_section_area
        self.volume = lenght*cross_section_area
        self.weight = self.volume * density

    def calculate(self, displacements_vector):
        displacements_vector = [displacements_vector[i] for i in self.liberty_degree]
        # return deformation, tension, interna forces
        deformation = (1/self.lenght) * np.array([-self.cos, -self.sin, self.cos, self.sin]).dot(displacements_vector)
        tension = self.modulus_of_elasticity * deformation
        internal_forces = tension * self.cross_section_area
        return deformation, tension, internal_forces
    

class Structure:

    def __init__(self):
        self.solved = False

    def load(self, fille_name):
        data = importa(fille_name)
        self.nodes_number = data[0]
        self.nodes_matrix = data[1]
        self.members_number = data[2]
        self.incidence_matrix = data[3]
        self.loads_number = data[4]
        self.loads_vector = data[5]
        self.restrictions_number = data[6]
        self.restrictions_vector = [int(x) for x in data[7]]
        self.colapse_conditions = data[8]
        self.solved = False

        self.nodes = []
        for node in range(self.nodes_number):
            self.nodes.append(Node(node, self.nodes_matrix[0][node], self.nodes_matrix[1][node]))
        
        self.members = []
        
        
        self.member_lenghts = np.zeros((self.members_number, 1))
        self.member_weights = np.zeros((self.members_number, 1))
        self.weight = np.zeros((1, 1))
        for i in range(self.members_number):
            member = Member(self.nodes[(int(self.incidence_matrix[i][0]))-1], self.nodes[(int(self.incidence_matrix[i][1]))-1], self.incidence_matrix[i][2], self.incidence_matrix[i][3], self.incidence_matrix[i][4])
            self.member_lenghts[i][0] = member.lenght
            self.member_weights[i][0] = member.weight
            self. weight[0,0] += member.weight
            self.members.append(member)

        self.members_number = len(self.members)

    def plot(self):
        plota(self.nodes_matrix, self.incidence_matrix)

    def solve(self):

        global_rigidity_matrix = np.zeros((self.nodes_number * 2, self.nodes_number * 2))
        for member in self.members:
            for line in range(4):
                for column in range(4):
                    global_rigidity_matrix[member.liberty_degree[line]][member.liberty_degree[column]] += member.rigidity_matrix[line][column]

        #Countour Conditions
        contour_loads_vector = np.delete(self.loads_vector, self.restrictions_vector, 0)
        contour_global_rigidity_matrix = np.delete(global_rigidity_matrix, self.restrictions_vector, 0)
        contour_global_rigidity_matrix = np.delete(contour_global_rigidity_matrix, self.restrictions_vector, 1)
        
        #Solve for U
        contour_displacements_vector = np.linalg.solve(contour_global_rigidity_matrix, contour_loads_vector)

        #Global U
        displacements_vector = np.zeros((self.nodes_number * 2, 1))
        index = np.arange(0, self.nodes_number*2, 1)
        index = np.delete(index, self.restrictions_vector, 0)
        for i in range(len(index)):
            displacements_vector[index[i]] = contour_displacements_vector[i]


        for i in range(self.nodes_number): # Deslocamentos multiplicados por 100 para visualização
            self.nodes_matrix[0][i] += displacements_vector[i*2] * 100
            self.nodes_matrix[1][i] += displacements_vector[1 + i*2] * 100
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

        displacements_vector = abs(displacements_vector)
        self.displacements_colapse = np.greater_equal(displacements_vector, self.colapse_conditions[0])

        tensions_vector = abs(tensions_vector)
        self.rupture_colapse = np.greater_equal(tensions_vector, self.colapse_conditions[1])

        percentage_deformation = abs(deformations_vector) / self.member_lenghts
        self.deformation_colapse = np.greater_equal(percentage_deformation, self.colapse_conditions[3])

        self.solved = True

    def write(self, fille_name):

        if self.solved:
            geraSaida(fille_name, [self.target_loads_vector, self.resultant_loads_vector], self.displacements_vector, self.displacements_colapse, self.deformations_vector, self.deformation_colapse, self.internal_forces_vector, self.tensions_vector, self.rupture_colapse, self.member_weights, self. weight, self.member_lenghts)
        else:
            print("Solve first")

structure = Structure()

structure.load("entrada_ponte") 
structure.plot()
structure.solve()
structure.plot()
structure.write("entrada_ponte")