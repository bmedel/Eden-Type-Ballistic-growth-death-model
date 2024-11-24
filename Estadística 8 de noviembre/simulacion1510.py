import funciones_new1510 as fn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import networkx as nx
import pickle

# Parámetros del sistema

parametros = [(1, 0), (0.9,0.1), (0.7, 0.3)]

estadistica_largos_param = []

for r_a, r_s in parametros:

    print([r_a,r_s])
    
    # Tamaños del sistema
    
    L0 = 6

    L = [L0, 2*L0, 4*L0]

    estadistica_largos_Li = []

    for m, size in enumerate(L):

        print("Tamaño = " + str(size))
        estadistica_largos = []
        # ------------------------------
        if m <= 1:
            realizaciones = 60
        else:
            realizaciones = 25
        for i in range(realizaciones):
            print(i)
            ind_activo = [0]
            ang = np.pi/3

            opt1 = np.array([np.cos(0*ang),np.sin(0*ang)])
            opt2 = np.array([np.cos(3*ang),np.sin(3*ang)])
            opt3 = np.array([np.cos(2*ang),np.sin(2*ang)])
            opt4 = np.array([np.cos(4*ang),np.sin(4*ang)])

            std = 1
            #------------------------------------------------

            positions_list = np.array([0*opt1])
            edgelist = [[0]]
            N_part = 1
            t_vida = np.ones(len(positions_list))
            t_crit = 8
            t_final = 20
        
            # Tamaño sistema
            bordes = (size, size)

            # Tamaño máximo del sistema

            N_max = size**2 + (3/2) * size + 1

            # Inicialmente 

            N_a = N_part
            N_s = 0

            while N_part/N_max<0.9:
                
                
                #----------------SELECCION DE PROCESO---------------------------------

                #------ Gillespie------------
                alpha_1 = r_a*N_a
                alpha_2 = r_s*N_s 
                alpha = alpha_1 + alpha_2

                p = alpha_1/alpha
                
                # Seleccion de numeros aleatorios
                
                u1 = np.random.rand()
                u2 = np.random.rand()

                # Seleccion de tiempo (gillespie)
                t_jump = 1/alpha * np.log(1-u2)
                
                if u1 < p:

                    #-----------------PROCESO DE CRECIMIENTO------------------------------
                    ##################----------------------###########################
                    a = fn.move(positions_list, ind_activo, edgelist, std, bordes)
                    positions_list = a
                    b = fn.check_vecinos(positions_list, edgelist, std, bordes)
                    ind_activo = b
                    N_part = N_part + 1
                    t_vida_el = 0
                    t_vida = np.append(t_vida, t_vida_el)
                    t_vida = t_vida - t_jump
                
                elif u1<1:
                    #-------------------PROCESO DE MUERTE---------------------------------
                    ##################----------------------############################
                    a = fn.kill(propensas, positions_list, t_vida, edgelist)
                    positions_list = a
                    b = fn.check_vecinos(positions_list, edgelist, std, bordes)
                    ind_activo = b
                    
                    N_part = N_part - 1
                    t_vida = t_vida - t_jump
            
                
                
                G = nx.DiGraph()
                G.add_edges_from([x for x in edgelist[1:] if isinstance(x, list)])
                
                puntas = {n for n, d in G.degree() if d == 1}
                

                propensas = []
                for terminal in puntas:
                    if t_vida[terminal] <= t_crit:
                        propensas.append(terminal)

                propensas = [x for x in propensas if x not in ind_activo]
                N_a = len(ind_activo)
                N_s = len(propensas)
            

            ######################### MEDIDAS DE LOS LARGOS DE LOS SEGMENTOS ################################### 

            largo_segmentos = []

            G = nx.DiGraph()
            G.add_edges_from([x for x in edgelist[1:] if isinstance(x, list)])
            ramificaciones = {n for n, d in G.degree() if d > 2}    
            puntas = {n for n, d in G.degree() if d == 1}
            nodos_intermedios = {n for n, d in G.degree() if d == 2}
            
            segmentos = []
            visitados = []
            for ramificacion in ramificaciones:
                for vecino in G.neighbors(ramificacion):
                    cont = 1
                    ruta = [ramificacion]  # Iniciamos la ruta con la ramificación
                    while vecino in nodos_intermedios:  # Exploramos hasta una punta o nueva ramificación
                        ruta.append(vecino)
                        for vecinos_aledaños in G.neighbors(vecino):
                            cont += 1
                            if vecinos_aledaños != vecino:
                                visitados.append(vecino)
                                vecino = vecinos_aledaños
                                
                    ruta.append(vecino)  # Agregamos el vecino final (punta o ramificación)
                    if vecino not in visitados:
                        largo_segmentos.append(cont)
                        segmentos.append(ruta)  # Guardamos toda la ruta
            
            estadistica_largos.append(largo_segmentos)

        estadistica_largos_Li.append(estadistica_largos)

    estadistica_largos_param.append(estadistica_largos_Li)


# # Nombre del archivo CSV
# filename = 'estadistica1311.pkl'

# # Abrir el archivo en modo escritura
# with open(filename, 'wb') as fp:
#     pickle.dump(estadistica_largos_param, fp)
    
    
        


