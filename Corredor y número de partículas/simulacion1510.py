import funciones_new1510 as fn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import networkx as nx
import copy

# import cProfile, pstats, io    # Este es un snippet de código para hacer profiling, ver cuanto y qué demora.
# profiler = cProfile.Profile()
# profiler.enable()

# -------------------------------


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
nodos_activos = []
N_part = 1
t_vida = np.ones(len(positions_list))
t_crit = 8
t_final = 300

# Tamaño del corredor
bordes = (22,3)

r_a = 0.9
r_s = 0.1

# Inicialmente 

N_a = N_part
N_s = 0

# Listas de estado

positions_list_estado = []
edgelist_estado = []
ind_activo_estado = []

while t_vida[0]<t_final:
    
    
    #----------------SELECCION DE PROCESO---------------------------------

    #----------------------Gillespie--------------------------------------
    alpha_1 = r_a*N_a
    alpha_2 = r_s*N_s 
    alpha = alpha_1 + alpha_2

    p = alpha_1/alpha
    
    # Seleccion de numeros aleatorios
    
    u1 = np.random.rand()
    u2 = np.random.rand()

    # Seleccion de tiempo (gillespie)
    t_jump = -1
    
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


    # Guardado de estado
    positions_list_estado.append(positions_list)
    edgelist_estado.append(copy.deepcopy(edgelist))
    ind_activo_estado.append(ind_activo)


        


#--------------------------------ANIMACIONES-----------------------------
#Crear la figura para la animación


from matplotlib.patches import Rectangle
import matplotlib.cm as cm  # Colormaps exóticos

fig, ax = plt.subplots()
ax.set_xlim(-bordes[0], bordes[0])
ax.set_ylim(-bordes[0], bordes[0])

colormap = cm.plasma

def update(i):
    ax.clear()
    ax.set_xlim(-bordes[0], bordes[0])
    ax.set_ylim(-bordes[0], bordes[0])
    positions_list_i = positions_list_estado[i]
    edgelist_i = edgelist_estado[i]
    ind_activo_i = ind_activo_estado[i]

    ax.scatter(positions_list_i[0,0], positions_list_i[0,1],color='red', s=20)
    #ax.scatter(positions_list_i[1:,0], positions_list_i[1:,1],color='black', s=5)
    ax.scatter(positions_list_i[ind_activo_i, 0], positions_list_i[ind_activo_i, 1], color='red', s=5)
    for j in range(1, len(edgelist_i)):
        if isinstance(edgelist_i[j], list):
            idxanterior = edgelist_i[j][0]
            idxsiguiente = edgelist_i[j][1]

            x_values = [positions_list_i[idxanterior][0], positions_list_i[idxsiguiente][0]]
            y_values = [positions_list_i[idxanterior][1], positions_list_i[idxsiguiente][1]]
            ax.plot(x_values, y_values, 'k-',color='black')

    # Crear un rectángulo centrado en (0, 0) con dimensiones proporcionadas por 'bordes'
    color_exotico = 'blue'  # Color basado en el frame
    rect = Rectangle(
        (-bordes[0], -bordes[1]),  # Esquina inferior izquierda
        2 * bordes[0], 2 * bordes[1],  # Ancho y alto (simétricos)
        linewidth=2, edgecolor=color_exotico, facecolor='none'
    )
    ax.add_patch(rect)

#Crear la animación
tiempo = 8000
intervalo = tiempo/len(positions_list_estado)
ani = FuncAnimation(fig, update, frames=len(positions_list_estado), interval=intervalo, repeat=True)

#Guardar la animación como un archivo .mp4 o mostrarla
#ani.save('corredorconmuerte1811(3).mp4', writer='ffmpeg')  # Para guardar la animación
plt.show()  # Para mostrar la animación
