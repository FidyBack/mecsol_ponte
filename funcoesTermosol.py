# -*- coding: utf-8 -*-
"""
A funcao 'plota' produz um gráfico da estrutura definida pela matriz de nos N 
e pela incidencia Inc.

Sugestao de uso:

from funcoesTermosol import plota
plota(N,Inc)
-------------------------------------------------------------------------------
A funcao 'importa' retorna o numero de nos [nn], a matriz dos nos [N], o numero
de membros [nm], a matriz de incidencia [Inc], o numero de cargas [nc], o vetor
carregamento [F], o numero de restricoes [nr] e o vetor de restricoes [R] 
contidos no arquivo de entrada.

Sugestao de uso:
    
from funcoesTermosol import importa
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xlsx')
-------------------------------------------------------------------------------
A funcao 'geraSaida' cria um arquivo nome.txt contendo as reacoes de apoio Ft, 
deslocamentos Ut, deformacoes Epsi, forcas Fi e tensoes Ti internas. 
As entradas devem ser vetores coluna.

Sugestao de uso:
    
from funcoesTermosol import geraSaida
geraSaida(nome,Ft,Ut,Epsi,Fi,Ti)
-------------------------------------------------------------------------------

"""
def plota(N,Inc):
    # Numero de membros
    nm = len(Inc[:,0])
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt

#    plt.show()
    fig = plt.figure()
    # Passa por todos os membros
    for i in range(nm):
        
        # encontra no inicial [n1] e final [n2] 
        n1 = int(Inc[i,0])
        n2 = int(Inc[i,1])        

        plt.plot([N[0,n1-1],N[0,n2-1]],[N[1,n1-1],N[1,n2-1]],color='r',linewidth=3)


    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.grid(True)
    plt.axis('equal')
    plt.show()
    
def importa(entradaNome):

    entradaNome = "entradas/" + entradaNome + ".xlsx"
    
    import numpy as np
    import xlrd
    
    arquivo = xlrd.open_workbook(entradaNome)
    
    ################################################## Ler os nos
    nos = arquivo.sheet_by_name('Nos')
    
    # Numero de nos
    nn = int(nos.cell(1,3).value)
                 
    # Matriz dos nós
    N = np.zeros((2,nn))
    
    for c in range(nn):
        N[0,c] = nos.cell(c+1,0).value
        N[1,c] = nos.cell(c+1,1).value
    
    ################################################## Ler a incidencia
    incid = arquivo.sheet_by_name('Incidencia')
    
    # Numero de membros
    nm = int(incid.cell(1,5).value)
                 
    # Matriz de incidencia
    Inc = np.zeros((nm,5))
    
    for c in range(nm):
        Inc[c,0] = int(incid.cell(c+1,0).value)
        Inc[c,1] = int(incid.cell(c+1,1).value)
        Inc[c,2] = incid.cell(c+1,2).value
        Inc[c,3] = incid.cell(c+1,3).value
        Inc[c,4] = incid.cell(c+1,4).value
    
    ################################################## Ler as cargas
    carg = arquivo.sheet_by_name('Carregamento')
    
    # Numero de cargas
    nc = int(carg.cell(1,4).value)
                 
    # Vetor carregamento
    F = np.zeros((nn*2,1))
    
    for c in range(nc):
        no = carg.cell(c+1,0).value
        xouy = carg.cell(c+1,1).value
        GDL = int(no*2-(2-xouy)) 
        F[GDL-1,0] = carg.cell(c+1,2).value
         
    ################################################## Ler restricoes
    restr = arquivo.sheet_by_name('Restricao')
    
    # Numero de restricoes
    nr = int(restr.cell(1,3).value)
                 
    # Vetor com os graus de liberdade restritos
    R = np.zeros((nr,1))
    
    for c in range(nr):
        no = restr.cell(c+1,0).value
        xouy = restr.cell(c+1,1).value
        GDL = no*2-(2-xouy) 
        R[c,0] = GDL-1

    # Vetor com as condições de colapso
    C = np.zeros((4,1))
    
    C[0][0] = restr.cell(1,2).value
    C[1][0] = restr.cell(1,4).value
    C[2][0] = restr.cell(1,5).value
    C[3][0] = restr.cell(1,6).value


    return nn,N,nm,Inc,nc,F,nr,R,C

def geraSaida(nome,TFt,Ut,CUt,Epsi,CEpsi,Fi,Ti,CTi,Pb,P,Tb):
    from xlutils.copy import copy # http://pypi.python.org/pypi/xlutils
    from xlrd import open_workbook # http://pypi.python.org/pypi/xlrd

    Ft = TFt[0]
    nome0 = "saidas/" + nome + '.txt'
    f = open(nome0,"w+")
    f.write('Reacoes de apoio [N]\n')
    f.write(str(Ft))
    f.write('\n\nDeslocamentos [m]\n')
    f.write(str(Ut))
    f.write('\n\nColapso por deslocamento [m]\n')
    f.write(str(CUt))
    f.write('\n\nDeformacoes []\n')
    f.write(str(Epsi))
    f.write('\n\nColapso por deformação [m]\n')
    f.write(str(CEpsi))
    f.write('\n\nForcas internas [N]\n')
    f.write(str(Fi))
    f.write('\n\nTensoes internas [Pa]\n')
    f.write(str(Ti))
    f.write('\n\nColapso por ruptura [m]\n')
    f.write(str(CTi))
    f.write('\n\nPeso das barras [kg]\n')
    f.write(str(Pb))
    f.write('\n\nPeso total [kg]\n')
    f.write(str(P))
    f.write('\n\nTamanho das barras [m]\n')
    f.write(str(Tb))
    f.close()
    
    Ft = TFt[1]
    dist = 10*" "
    nome1 = "saidas/" + nome + "_format" + '.txt'
    f = open(nome1,"w+")
    f.write('Reacoes de apoio [N] e nDeslocamentos [m]\n')
    n = int(len(Ft) / 2)
    for i in range(n):
        x = "{:.2f}".format(Ft[0 + i*2][0])
        y = "{:.2f}".format(Ft[1 + i*2][0])
        f.write("Nó {:02d} x: {:>9} [N] y: {:>9} [N]{}x: {:+.2e} [m] y: {:+.2e} [m]\n".format(i+1, x, y, dist, Ut[0 + i*2][0], Ut[1 + i*2][0]))

    dist = 3*" "
    f.write('\nDeformacoes [], Forcas internas [kN], Tensoes internas [GPa], Peso das barras [g] e Tamanho das barras [mm]\n')
    for u in range(len(Epsi)):
        a = "{:.2f}".format(Fi[u][0]/1000)
        b = "{:.2f}".format(Ti[u][0]/1000000)
        c = "{:.2f}".format(Pb[u][0]*1000)
        d = "{:.2f}".format(Tb[u][0]*1000)
        f.write("Elemento {:02d}: {:+.2e} []{}{:>9} [kN]{}{:>9} [GPa]{}{:>9} [g]{}{:>9} [mm]\n".format(u+1, Epsi[u][0], dist, a, dist, b, dist, c, dist, d))

    f.write('\nPeso total [g]\n')
    e = "{:.2f}".format(P[0][0]*1000)
    f.write("Peso: {:>7}\n".format(e))

    f.close()

    entradaNome = "entradas/" + nome + ".xlsx"
    original = open_workbook(entradaNome)
    copia = copy(original) 
    incid_sheet = copia.get_sheet('Incidencia') 

    minimo = abs(min(Ti))
    maximo = max(Ti)
    if(minimo > maximo):
        maximo = minimo

    for i in range(len(Ti)):
        if Ti[i][0] > 0:
            incid_sheet.write(i+1, 7, float(Ti[i][0]/maximo))
        else:
            incid_sheet.write(i+1, 7, float(Ti[i][0]/-maximo))

    incid_sheet.write(1, 8, float(P[0][0]))

    

    copia.save("entradas/" + nome + ".xlsx")

