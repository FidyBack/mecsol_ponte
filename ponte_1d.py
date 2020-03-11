import numpy as np

# Variáveis
elem = 2
L = 2
A = 0.02
P = 50*(10**3)
E = 200*(10**9)

# Calcula os vetores de deslocamento nodal de acordo com a quatidade de elementos passados
def vet_des_nod(elem, L, A, P, E):
    # Matriz de Rigidez
    for i in range(elem):
        Mf = np.zeros((elem+1, elem+1))
    for j, element in enumerate(Mf):
        for k1, element in enumerate(Mf):
            if (k1 == j+1 or k1 == j-1):
                Mf[j][k1] = -1
            if (k1 == j and k1 !=0 and k1 != elem):
                Mf[j][k1] = 2
            if (k1 == j == 0 or k1 == j == elem):
                Mf[j][k1] = 1
    k = ((E*A)/(L/elem))*Mf

    # Vetor de Forças
    F = np.zeros((elem+1, 1))
    for j, element in enumerate(F):
        if (j == elem):
            F[j] = P

    # Condições de Contorno
    Fcc = np.delete(F, 0, 0 )
    Klin = np.delete(k, 0, 0)
    Kcc = np.delete(Klin, 0, 1)

    # Inversão e Multiplicação
    U = np.linalg.solve(Kcc, Fcc)
    return(U, k)

# Calcula a reação de apoio no nó 1
def reac_apoio(vet_nod):
    mat_vet = np.vstack (([0], vet_nod[0])) 
    mat_k = vet_nod[1]
    mat_multip = np.matmul(mat_k, mat_vet)
    return mat_multip

# Calcula a deformação longitudinal específica em cada elemento
def def_esp(vet_nod, L, elem):
    mat_vet = np.vstack (([0], vet_nod[0]))
    mat_def = np.zeros((elem, 1))
    for i in range(len(mat_vet)-1):
        def1 = mat_vet[i+1]-mat_vet[i]/(L/elem)
        mat_def[i] = def1
    return mat_def

# Calcula a Tensão em Cada elemento
def tens_elem(def_esp, E):
    return def_esp*E

vet_nod = vet_des_nod(elem, L, A, P, E)
forc_apoio = reac_apoio(vet_nod)
deform = def_esp(vet_nod, L, elem)
tens = tens_elem(deform, E)

# Prints
print('-----------------------')
print('Vetor de Deslocamento Nodal [m]')
print("{}".format(vet_nod[0]))

print('-----------------------')
print("Reação de apoio no nó 1 [N]")
print("{}".format(forc_apoio[0]))

print('-----------------------')
print("Deformação longitudinal em cada elemento")
print("{}\n".format(deform))

print('-----------------------')
print("Tensão em cada elemento")
print("{}\n".format(tens))