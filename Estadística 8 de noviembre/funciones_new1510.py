import numpy as np
import matplotlib.pyplot as plt
import random as rand
from matplotlib.animation import FuncAnimation
from scipy.optimize import curve_fit
import scipy.stats as sta
import networkx as nx


# Este es el archivo con las funciones buenas (18 de noviembre).

def direccion_priv(pos_part_ant, part_nueva,std, options):
    """Dado un par de partículas dirigidas, podemos encontrar la dirección privilegiada, es decir, 
    el índice asociado a la dirección que tomó la partícula anterior. Esta función devuelve el arreglo normalizado
    de las probabilidades de las opciones considerando una distribución en especifico."""
   
    for index, lugar in enumerate(options):
        if np.allclose(pos_part_ant + lugar, part_nueva,rtol=1e-05, atol=1e-08, equal_nan=False):
            index_priv = index
            break
            
        else:
            index_priv = None
    
  
    if index_priv + 3 >= len(options):
      options = np.delete(options, index_priv - 3,axis=0)
      index_priv = index_priv - 1
    else:
      options = np.delete(options, index_priv + 3,axis=0)
    
    # Tenemos bien la lista de opciones.
    indices = np.arange(len(options))
    
      
    distancias = np.minimum(np.abs(indices - index_priv), len(options) - np.abs(indices - index_priv))
    array_prob = sta.norm.pdf(distancias, loc = 0, scale=std)

    array_prob = array_prob/array_prob.sum()

       
    return array_prob, options

def move(positions_list, ind_activo,edgelist, std, bordes):
      """Función que toma una lista de posiciones con un elemento inicial, chequea el entorno y elige aleatoriamente entre espacios libres. Luego, añade el índice
      de la partícula (cada índice sería el largo de la lista posiciones) a una lista de índices"""
      # Primero realizamos un test sobre el entorno a la hora de moverse.
      # Primero arreglamos el error de seleccionar aleatoriamente entre un elemento (raro)
      ang = np.pi/3
      options=np.array([(np.cos(0*ang),np.sin(0*ang)),(np.cos(ang),np.sin(ang)),(np.cos(2*ang),np.sin(2*ang)),
                        (np.cos(3*ang),np.sin(3*ang)), (np.cos(4*ang),np.sin(4*ang)), 
                        (np.cos(5*ang),np.sin(5*ang))])
      
      
      indice = rand.choice(ind_activo)
      # Construiremos un árbol genealógica de modo que cada posición tenga contacto con sus ancestros
      if len(edgelist) == 1:
         ancestros = [0]
      else:
        ancestros = [indice]

      edgelist.append(ancestros + [len(positions_list)])
      #--------------------------------------------------------------------------
      # IMPLEMENTACION DEL MODELO BALISTICO
      
      # Asociamos a la posición
      position_actual = positions_list[indice]
      
      # Haremos una serie de tests...
      
      if not len(edgelist[indice]) == 2:
        
        # Exploramos las opciones y descartamos aquellos lugares donde haya otra partícula
        indexrip =[]
        for i, option in enumerate(options):
          lugar_hipotetico= position_actual + option
          test = np.all(np.isclose(lugar_hipotetico, positions_list, rtol=1e-05, atol=1e-08, equal_nan=False), axis=1)
          
          if np.any(test):
            indexrip.append(i)
          else:
            pass
      # Finalmente eliminamos espacios que ya están ocupados
      
        options = np.delete(options, indexrip,axis=0)
        j = np.random.randint(len(options))
        nueva_pos = np.array([(position_actual + options[j])]) # Se elige una de estas posiciones aleatoriamente
      
      else:
        
        # Vemos el ancestro directo de la partícula
        ind_position_pasada = edgelist[indice][0]
        
        position_pasada = positions_list[ind_position_pasada]

        # Aquí tenemos que añadir la dirección:
        # Esto se puede realizar determinando el lugar elegido anteriormente, cómo hacemos esto?
        # podemos encontrar la dirección en la que está la partícula anterior.

        array_prob, options_new = direccion_priv(position_pasada, position_actual, std, options)

        # Notar bien que index_rip es el nuevo indexrip con la posición ocupada trivial considerada.
        # Ahora, lo que queremos hacer es acomodar la distribución tal que esté centrada en el indice privilegiado.

        # Exploramos las opciones y descartamos aquellos lugares donde haya otra partícula
    
        indexrip =[]
        for i, option in enumerate(options_new):
          lugar_hipotetico= position_actual + option
          test = np.all(np.isclose(lugar_hipotetico, positions_list, rtol=1e-05, atol=1e-08, equal_nan=False), axis=1)
          test1 = np.all(np.abs(lugar_hipotetico) <= np.array([bordes]))
          if np.any(test) or not test1:
            indexrip.append(i)
        # Finalmente eliminamos espacios que ya están ocupados
        options_new = np.delete(options_new, indexrip,axis=0)
        array_prob = np.delete(array_prob, indexrip, axis=0)
        array_prob_norm = array_prob/np.sum(array_prob)
        
        # Finalmente, escogemos

        j = np.random.choice(len(options_new), 1, p=array_prob_norm)
        nueva_pos = np.array([(position_actual + options_new[j[0]])])
        
        # Y se las damos a las posiciones
      
      positions_list = np.append(positions_list, nueva_pos, axis=0) # Se agrega a la lista de posiciones
      return positions_list

def check_vecinos(positions_list, edgelist, std, bordes):
    """Toma la lista de posiciones y chequea el entorno de cada partícula individual.
    Si hay espacios libres, no hace nada. Si no hay espacios libres, elimina a aquella partícula de las partículas activas."""
    
    ind_activo = []  # Reiniciamos la lista de partículas activas
    ang = np.pi / 3  # Ángulo de 60 grados
    options = np.array([(np.cos(0 * ang), np.sin(0 * ang)),
                        (np.cos(1 * ang), np.sin(1 * ang)),
                        (np.cos(2 * ang), np.sin(2 * ang)),
                        (np.cos(3 * ang), np.sin(3 * ang)),
                        (np.cos(4 * ang), np.sin(4 * ang)),
                        (np.cos(5 * ang), np.sin(5 * ang))])

    for i, j in enumerate(positions_list): # Trabajaremos sobre cada partícula

      if np.any(np.isnan(j)): # Testea si no es una partícula muerta (positions_list[muerta] = NaN)
         pass
      
      else:
        if not len(edgelist[i]) == 2: # Si la partícula está en el inicio, no tendrá una dirección privilegiada (ergo equiprob)

          test = np.all([
              np.any(np.all(np.isclose(j + k, positions_list, rtol=1e-5, atol=1e-6, equal_nan=False), axis=1))
              for k in options]) # Testea si hay algun espacio en alguna dirección que no esté ocupado por alguna partícula en position_list.
          test2 = np.any([np.all(np.abs(j+k) <= np.array([bordes]))] for k in options) # Revisa si hay alguna dirección que se mantenga dentro de la caja (bordes)
          if not test and test2: # Si se cumplen estas condiciones, es partícula activa.
              ind_activo.append(i)
        else:
           ind_position_pasada = edgelist[i][0] # revisamos la partícula de la cual viene la partícula activa
        
           position_pasada = positions_list[ind_position_pasada]
           position_actual = positions_list[i]

            # Aquí tenemos que añadir la dirección:
            # Esto se puede realizar determinando el lugar elegido anteriormente, cómo hacemos esto?
            # podemos encontrar la dirección en la que está la partícula anterior.

           array_prob, options_new = direccion_priv(position_pasada, position_actual, std, options) # option_new devuelve el arreglo exceptuando la opción donde se encuentra la partícula pasada

           for k in options_new:
            test = np.any(np.all(np.isclose(j + k, positions_list, rtol=1e-5, atol=1e-6, equal_nan=False), axis=1)) # Mismos test que antes, notar que la misma opción debe cumplir ambas condiciones
            test1 = np.all(np.abs(j+k) <= np.array([bordes]))
            if not test and test1:
                ind_activo.append(i)
                break # Si alguna dirección cumple con esta condición, es suficiente.
    
    return ind_activo



def kill(propensas, position_list, life_times, edgelist):
  """Elimina del sistema la partícula seleccionada"""
  muerta=rand.choice(propensas)
  position_list[muerta]=np.nan
  edgelist[muerta] = None
  life_times[muerta] = np.inf
  return position_list
