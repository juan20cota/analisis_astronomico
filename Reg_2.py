# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 15:08:51 2022

@author: Juan Carlos De Borbón Cota 

ANALISIS DE DATOS DE REGIONES DE FORMACIÓN ESTELAR

En este programa usamos la paquetería de Astrodendro,
el cual es un progrma orientado al analizis de objetos 
astronómicos como la evolución estelar, morfología galáctica
y regiones de formación estelar, este último siendo nuestro interés
de estudio para el analisis de las propiedades espectrales de los
datos obtenidos del catálogo de MaNGA del SDSS.


----------------------------------------------------------------------
Para mas información acerca del programa astrodendro favor de revisar
la documentación en su página oficial: https://dendrograms.readthedocs.io/en/stable/index.html
Ahí encontraran manuales para la instalación y función del programa astrodendro.

"""
# Importación de las bibliotecas necesarias
import os
import glob
import matplotlib
from astrodendro.analysis import ScalarStatistic
import pandas as pd
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astrodendro import Dendrogram, pp_catalog
from astrodendro.analysis import PPStatistic
from astrodendro.pruning import min_peak
from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import LinearLocator
from matplotlib.path import Path
import matplotlib.cm as cm
from matplotlib import colormaps
import numpy as np
# Función para cargar archivos con extensión '.cube.fits' de una carpeta específica
def cargar_archivos(carpeta):
    archivos_fits = []  # Lista donde se almacenarán los nombres de los archivos FITS
    
    # Cambiar el directorio de trabajo a la carpeta indicada
    os.chdir(carpeta)
    
    # Buscar todos los archivos en la carpeta que terminen con la extensión '.cube.fits'
    archivos = glob.glob("*.cube.fits")
    
    return archivos  # Retorna la lista de archivos encontrados

# Función para cargar y procesar los archivos FITS
def cargar_fits(carpeta):
    # Llamamos a la función cargar_archivos para obtener los archivos .cube.fits
    archivos_cubefits = cargar_archivos(carpeta)
    
    datos_fits = []  # Lista donde se almacenarán los datos procesados de los archivos FITS
    
    # Iteramos sobre cada archivo FITS encontrado
    for archivo in archivos_cubefits:
        # Abrir el archivo FITS usando astropy
        img = fits.open(archivo)
        
        # Extraer los datos de la extensión "FLUX_ELINES" del archivo FITS
        lines = img["FLUX_ELINES"].data
        
        # Seleccionar una sección específica de los datos (por ejemplo, la línea H-alpha)
        h_alpha = lines[45, :, :]
        
        # Agregar los datos de H_alpha a la lista de datos
        datos_fits.append(h_alpha)
        
        # Cerrar el archivo FITS después de haberlo procesado para liberar recursos
        img.close()
    
    return datos_fits  # Retorna la lista con los datos procesados de todos los archivos

# Llamar a la función para cargar los archivos FITS y obtener sus datos
carpeta = r'C:\Users\juan2\OneDrive\Escritorio\Juan\servicio\astrdendro'
datos = cargar_fits(carpeta)



# =============================================================================
# CÓDIGO PARA EL ANÁLISIS CON ASTRODENDRO
# =============================================================================

# Función que realiza el análisis utilizando la librería astrodendro
def dendro_data(datos):

   
    gal = cargar_archivos(carpeta)
    pixel_values= []    

    

    # Iterar sobre cada conjunto de datos cargados
    for datos in datos:
        # Calcular el dendrograma usando los datos y configuraciones específicas
        d = Dendrogram.compute(datos, min_value=0.0, min_delta=0.0, verbose=True)
        
        # Crear el objeto plotter para visualizar el dendrograma
        p = d.plotter()
        
        # Obtener las hojas (leaves) del dendrograma
        l = d.leaves
    
        # Lista para almacenar las alturas de las hojas
        maxv = []
        
        
        

        # Iterar sobre cada hoja del dendrograma
        for leaf in d.leaves:  
            # Agregar la altura de cada hoja a la lista maxv
            maxv.append(leaf.height)
        
        # Calcular el valor máximo de las alturas y dividirlo por 2
        MAX = (max(maxv)) / 3
        
        # Lista para almacenar los valores mayores que MAX
        valores_mayores = []
    
        # Iterar sobre cada valor en la lista maxv
        for valor in maxv:
            if valor > MAX:
                valores_mayores.append(valor)  # Guardar el valor si es mayor que MAX
    
        # Usar la función min_peak para obtener un umbral independiente
        custom_independent = min_peak(min(valores_mayores))
    
        # Volver a calcular el dendrograma utilizando el umbral personalizado
        d = Dendrogram.compute(datos, is_independent=custom_independent, verbose=True)
        
        # Crear el objeto plotter para la nueva visualización
        p = d.plotter()
    
        # Crear una figura de Matplotlib para mostrar los resultados
        fig = plt.figure()
    
        # Crear el primer subplot (lado derecho) para el árbol del dendrograma
        ax = fig.add_subplot(122)
        p.plot_tree(ax, color='black')  # Dibujar el árbol del dendrograma en negro
    
        # Crear el segundo subplot (lado izquierdo) para la imagen de los datos
        ax = fig.add_subplot(121)
        ax.imshow(datos, cmap='gray')  # Mostrar la imagen de los datos en escala de grises
    
        Reg = len(d.leaves)
    # for galaxias in gal:
        
    #     print(f"Número de regiones en la imagen {galaxias}: {Reg}")
         
        # Crear una paleta de colores

        cmap =colormaps['tab10'].resampled(Reg)  # Reg es el número de colores requeridos
        
    
        
        # Iterar sobre cada hoja del dendrograma
        for leaf_idx, leaf in enumerate(d.leaves, start=1):
           
            print(f"Procesando región {leaf_idx}")
    
            # Asignar un color único basado en el índice
            # Obtener el color del colormap normalizado
            rgba_color = cmap(leaf_idx - 1)  # Usar un índice entero para colores discretos
            color_hex = matplotlib.colors.to_hex(rgba_color)  # Convertir a formato hexadecimal
            
            print(f"Color para la región {leaf_idx}: {color_hex}")  # Mostrará algo como '#440154'
        
            # Dibujar el contorno de la estructura de cada hoja en la imagen
            mask, pv = p.plot_contour(ax, structure=leaf, color=color_hex)
            
            pixel_values.append(pv)
            #print(pixel_values)
            
            # Encontrar el centroide aproximado de la región
            ny, nx = mask.shape
            y_coords, x_coords = np.where(mask > 0)
            centroid_x = np.mean(x_coords)
            centroid_y = np.mean(y_coords)
        
            # Etiquetar la región con su índice
            ax.text(
                centroid_x, centroid_y, str(leaf_idx),
                color='black', fontsize=10, ha='center', va='center',
                bbox=dict(facecolor='lightgray', edgecolor='none', alpha=0.6, pad=0)  # Fondo para visibilidad
            )
        
        print("Proceso completado.")
        print("===================================================================")    
            
    return  pixel_values 





"""            
            # print(contours.collections[0].get_paths())    #ESTA JODIDA COSA DEBERÍA DE FUNCIONAAAAAR!
            # # Verificar si se generaron contornos
            # if contours is not None and len(contours.collections) > 0:
            #     paths = contours.collections[0].get_paths()
            #     print(f"Se encontraron {len(paths)} paths para la hoja {leaf}.")
            # else:
            #     print(f"No se generaron contornos para la hoja {leaf}.")
"""            
"""           
            ESTA BASURA NO SIRVIÓ, ME DA LOS VALORES DE LOS PIXELES PERO EN UN ARRAY DE 1D!, ASÍ QUE CHISTE
            MEJOR SIGO CON MI MÉTODO :(
                
            UPDATE: NO PUEDE SER!! MI MÉTODO TAMBIEN ME DA UN ARRAY DE 1D :( ESO QUIERE DECIR QUE TAL VEZ LA FUNCIÓN
                    DE ASTRODENDRO FUNCIONE, YA VEREMOS CUAL ESTÁ MEJOR, YO ESPERO QUE MI MÉTODO SEA MEJOR,
                    PORQUE LA MANERA EN LA QUE CALCULO EL ÁREA ES CON LA FÓRMULA DE POLÍGONOS Y
                    ASTRODENDRO LO HACE CONTANDO LOS PIXELES.

                    LA CUESTIÓN ES LA SIGUIENTE: LOS VALORES DE ASTRODENDRO SON MAYORES QUE 1,
                    MIENTRAS QUE EL METODO CON AX.CONTOUR ME DA LOS VALORES DEL PRIMER CONTORNO (INCLUYE 0).
                    
                    POR LO TANTO, AL MOMENTO DE CALCULAR EL ÁREA Y OBTENER UN VALOR PROMEDIO DE LOS PIXELES,
                    AMBOS MÉTODOS TENDRAN RESULTADOS DISTINTOS, LO CUAL ME PREOCUPA, PODRÍA HACER UN FILTRO 
                    PODRÍA HACER UN FILTRO QUE REMUEVA LOS PIXELES CON VALOR 0. CREO QUE ESA ES LA RESPUESTA                                                    
            
           UPDATE 2: VALIÓ ÑONGA, HICE EL FILTRO Y LOS VALORES MÍNIMOS Y MÁXIMOS DE LOS ARRAY SON IGUALES, FAK.  
                     YA VEREMOS SI LOS PROMEDIOS SON IGUALES, YO ESPERARÍA QUE NO.
    print("######################################################################")
    
    stat = PPStatistic(leaf)
    print(stat.area_exact,"area",stat.area_exact.value, "valor") CACAAAAAAAAAAAAA!  #UPDATE: NI TANTO, LO SUBESTIMÉ :(
    print(stat.stat.values)
    print("######################################################################") 
    
    
    stot = stat.stat.values 
    
    
    
    area_kk = stat.stat.count() * stot


            
            
"""
"""
            stat = PPStatistic(leaf)
            #stot = ScalarStatistic(leaf.stats.value, leaf.stat.indices)
            print(stat.area_exact,"area",stat.area_exact.value, "valor")    #CÓMO COÑO OBTENGO LOS VALORES DE LOS PIXELES?!!
            print(stat.x_cen,"x", stat.x_cen.value, "valor")
            print(stat.y_cen,"y")
            #print(stot.values, "intensidad según")
     """       
"""  
    PURA CACA DE ESTADÍSTICA     


   NAAAA NO CREO QUE ME SIRVA LO DE ARRIBA POR QUE NI SIQUIERA SÉ SI LOS DATOS
   QUE DAN SON DE LAS REGIONES ENCERRADAS EN ROJO PORQUE USAN EL COMANDO PLT.COUNTOUR PARA HACER LOS CONTORNOS,
   ASÍ QUE TENGO QUE BUSCAR LA MANERA EN QUE LA FUNCIÓN CONTOUR ME DÉ LOS VALORES Y 
   PARÁMETROS QUE USA PARA OBTENER LOS VALORES,
   DESPUÉS TENGO QUE APLICAR ESA MISMA REGIÓN A LAS DIFERENTES LÍNEAS, ASÍ QUE TENGO QUE BUSCAR LOS DATOS ESPACIALES,
   PARA INTRODUCIRLOS EN LA FUNCIÓN CONTOUR Y QUE ME DE LAS MISMAS REGIONES PARA LAS DIFERENTES LONGITUDES DE ONDA, 
   POR EL MOMENTO SOLO TENGO QUE HACER QUE FUNCIONE CON LA LÍNEA DE H-ALPHA, PERO DEBO TENER CUIDADO PARA QUE NO DELIMITE
   REGIONES NUEVAS SEGÚN LAS INTENSIDADES DE LAS DEMÁS LÍNEAS, ASÍ NO FUNCIONARÍA. 
   
   
"""            
      

"""
QUE NO SE TE OLVIDE USAR LOS DATOS DE LOS CENTROS, GUARDALOS EN UNA LISTA Y ESOS DATOS LOS USAREMOS PARA LAS DEMÁS LÍNEAS
SABRÁ DIOS CÓMO 

POR FAVOR DIOS AMPÁRAME :(

"""


# Llamar a la función dendro_data para realizar el análisis con los datos
val = dendro_data(datos)










# MISCELANEOUS





# # First image control

# fig, ha = plt .subplots()
# ha.imshow(h_alpha)


# # Second image of control
# fig, ha = plt.subplots()
# o = plt.plot(h_alpha)

# # Third image of control

# ax = plt.figure().add_subplot(111, projection='3d')

# x = np.arange(0, h_alpha[0].size)

# y = np.arange(0, h_alpha[1].size)

# X, Y = np.meshgrid(x, y)

# ax.plot_surface(X, Y, h_alpha, cmap=cm.viridis, edgecolor='g', lw=0.5, rstride=2, cstride=2, alpha=0.5)

# # Plot projections of the contours for each dimension.  By choosing offsets
# # that match the appropriate axes limits, the projected contours will sit on
# # the 'walls' of the graph.


# ax.contour(X, Y, h_alpha, zdir='z', offset=0, cmap='coolwarm', alpha=1)
# ax.contour(X, Y, h_alpha, zdir='x', offset=0, cmap='viridis')
# ax.contour(X, Y, h_alpha, zdir='y', offset=0, cmap='magma')



# #fourth image of control

# # Create a figure and plot the image
# fig, ax = plt.subplots()
# ax.imshow(h_alpha, cmap='gray')

# # Define the contour levels
# contour_levels = [1.8]

# contour = plt.contour(h_alpha, levels=[0.7], colors='red', linewidths=0.5)
# plt_label = plt.clabel(contour, inline=True, fontsize=8)
# plt.setp(plt_label, color='red')


# contour = plt.contour(h_alpha, colors='blue', linewidths=0.5)
# plt_label = plt.clabel(contour, inline=True, fontsize=8)
# plt.setp(plt_label, color='blue')

# #  Show the plot
# plt.show()








#gal 2
#1.324 = y              4.496 = x          n =  30%
# 
#


#gal 1
#0.785 = y              1.341 = x          n = 58.53%


#gal 1
#0.716                  0.716             
#0.783                  1.341             
#0.067 = y              0.625 = x          

#n = 10.72%



#gal 2
# base = 1.069
# y = 1.324 - base
# x = 4.496 - base
# n = 7.44%

# n = y*(1/(x/100))           



