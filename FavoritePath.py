import numpy as np
import math as ma
import matplotlib.pyplot as plt
from scipy.stats import loggamma

def matris(n,a):
    aux = []
    if a == 1:
        aux1 = np.random.normal(1 , 1, n*n)
    elif a == 2:
        aux1 = np.random.binomial(1 , 0.5, n*n)
    elif a == 3:
        aux1 = loggamma.rvs(1, loc=0, scale=1, size=(n*n), random_state=None)
    elif a == 4:
        aux1 = np.random.exponential(scale=10, size=(n*n))
    elif a == 5:
        aux1 = (np.random.pareto(10, (n*n))+1)*200
    for i in range(n):
        aux.append(aux1[i*n:(i+1)*n]) 
    return aux
    

#Encontraremos el peso de las trayectorias
#Recivimos un ambiente y regresamos el mayor de los pesos hasta ese punto
    
def ptrayect(m):
    amb = []        #guarda el ambiente con la suma de los pesos
    tij = []        #guarda las coordenadas de donte procede el punto anterior
    auxlen = len(m)
    for i in range(auxlen):
        auxamb = []
        auxtij = []
        for j in range(auxlen):
            if (i != 0 and j != 0):
                aux1 = amb[i-1][j] + m[i][j]
                aux2 = auxamb[j-1] + m[i][j]
                if aux1 > aux2:
                    auxamb.append(aux1)
                    auxtij.append([i-1,j])
                else:
                    auxamb.append(aux2)
                    auxtij.append([i,j-1])
            elif (i !=0 and j == 0):
                auxamb.append(amb[i-1][j] + m[i][j])
                auxtij.append([i-1,j])
            elif (i == 0 and j != 0 ):
                auxamb.append(auxamb[j-1] + m[i][j])
                auxtij.append([i,j-1])
            else:
                auxamb.append(m[i][j])
                auxtij.append([i,j])
        amb.append(auxamb)
        tij.append(auxtij)
    return(amb, tij)


# Creamos una función que nos regresa la masa y las coordenadas del punto final favorito de longitud n
def Loc_Masa_punto_final_fav(M,n):
    Loc_punto_final = []
    Masa_punto_final = []
    for i in range(n):
        aux_Loc_punto_final = []
        aux_Masa_punto_final = []
        if i != 0:
            for j in range(i):
                k = i-j
                aux_Loc_punto_final.append([j,k])
                aux_Masa_punto_final.append(M[j][k])
        else:
            aux_Loc_punto_final.append([0,0])
            aux_Masa_punto_final.append(M[0][0])
        Masa_punto_final.append(max(aux_Masa_punto_final))
        aux = aux_Masa_punto_final.index(max(aux_Masa_punto_final))
        Loc_punto_final.append(aux_Loc_punto_final[aux])     
    return(Masa_punto_final, Loc_punto_final)

#Funcion que simula m caminatas aleatorias
def CA(v,m):    #v ambiente, m numero de iteraciones
    C = []
    c = list(np.zeros((len(v),len(v))))
    c[0][0] = m
    for i in range(m):
        a = 0
        b = 0
        aux = np.random.binomial(1 , 0.5, len(v)-1)
        for j in aux:
            if j == 1:
                a += 1
            else:
                b += 1
            if [a,b] in C:
                c[a][b] += 1
            else:
                C.append([a,b])
                c[a][b] += 1
        if i % 1000 == 0:
            print(i)
    return(c)

#función que normaliza las diagonales
def Diag(v,m):      #v ambiente, c peso de la caminata
    c = CA(v,m)
    D = []
    for i in range(len(v)):
        aux = []
        for j in range(len(v)):
            aux.append(c[j][i-j]*v[j][i-j]/m)
        D.append(aux)
    return(D)



#funcion que resive coordenadas del punto favorito y las grafica.
def Grafica_punto_final(Lpf):
    auxx= []
    auxy= []
    for i in Lpf:
        auxx.append(i[0])
        auxy.append(i[1])
    plt.plot(auxx,auxy)
    plt.xlim(0, len(Lpf))
    plt.ylim(0, len(Lpf))
    plt.xlabel("X")
    plt.ylabel("Y")
    #plt.title("Localizacion Punto Favorito")
    plt.show()
    
def Grafica_Masa_punto_final(M,B):
    Z=np.zeros((len(M),len(M)))
    for i in range(len(M)):
        for j in range(len(M)):
            if i+j <= len(M):
                Z[i,j] = max(0,ma.floor((B*float(M[i][j]))/ma.log(10)))
            else:
                Z[i,j] = 0
    levels = [0, 0.0001 , 5, 10, 20, 40, 80, 100]
    plt.contourf(range(len(M)), range(len(M)), np.transpose(Z), levels=levels, extend='both', colors=('w','#ff0000', '#ff9900', '#999900', '#009999', '#0099ff', '#0000ff'))
    plt.xlabel("X")
    plt.ylabel("Y")
    #plt.title("Densidad Punto Favorito")
    plt.colorbar()
    plt.show()
    return()
    
#Funcion que itera los codigos anteriores para obtener la grafica del punto favorito
def IGraficasPF(n,a,B,ni,m):   #n tamaño del ambiente, a tipo de ambiente, B lista de temperaturas, ni numero de iteraciones, m numero de caminatas aleatorias
    for i in range(ni):
        v = matris(n,a)
        M, T = ptrayect(v)
        Mpf, Lpf = Loc_Masa_punto_final_fav(M,n)
        Grafica_punto_final(Lpf)
        print(Lpf[-1])
        D = Diag(v,m)
        for i in B:
            Grafica_Masa_punto_final(D,i)
    return()