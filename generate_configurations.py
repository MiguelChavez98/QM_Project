# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 13:52:21 2024

@author: migue
"""

import numpy as np
from scipy.spatial.transform import Rotation as R
from datetime import datetime


def xyz_a_matriz_numpy(ruta_archivo_xyz):
    # Leer el contenido del archivo
    with open(ruta_archivo_xyz, 'r') as archivo:
        lineas = archivo.readlines()

    # Omitir la primera línea (contiene el número total de átomos) y otras no necesarias
    # Asumimos que las coordenadas comienzan desde la tercera línea
    lineas_coordenadas = lineas[2:]
    
    # Preparar una lista para almacenar las coordenadas
    coordenadas = []
    
    for linea in lineas_coordenadas:
        if linea.strip():  # Verificar que la línea no esté vacía
            partes = linea.split()
            if len(partes) == 4:  # Asegurar que la línea tenga el formato esperado (Elemento x y z)
                # Convertir las partes numéricas (x, y, z) a flotantes y añadir a la lista
                coordenadas.append([float(partes[1]), float(partes[2]), float(partes[3])])

    # Convertir la lista de coordenadas a una matriz NumPy
    matriz_coordenadas = np.array(coordenadas)

    return matriz_coordenadas


def generate_random_points(distance_start, distance_end, num_spheres,
                           num_points_per_sphere, center=None):
    if center is None:
        center = np.array([0, 0, 0])

    all_points = []
    radius_values = np.linspace(distance_start, distance_end, num_spheres)

    for radius in radius_values:
        # Generación directa de puntos sobre la esfera
        points = np.random.normal(size=(num_points_per_sphere, 3))
        norms = np.linalg.norm(points, axis=1, keepdims=True)
        normalized_points = points / norms * radius
        translated_points = normalized_points + center
        all_points.append(translated_points)

    # Combinar y redondear en un solo paso
    all_points = np.around(np.vstack(all_points), decimals=8)

    return all_points




def random_points_sector(distance_start, distance_end, num_spheres,
                         num_points_per_sphere, center=None,
                         min_points_in_quadrant=None, max_points_in_quadrant=None):
    """
    Genera puntos aleatorios distribuidos uniformemente sobre la superficie de
    varias esferas concéntricas y verifica que cada capa cumpla con un número
    mínimo y máximo de puntos en el cuadrante positivo antes de continuar.
    """
    # Establecer valores predeterminados seguros para los parámetros
    if center is None:
        center = np.array([0, 0, 0])
    if min_points_in_quadrant is None:
        min_points_in_quadrant = 1  # Asumir al menos 1 punto como mínimo por defecto
    if max_points_in_quadrant is None:
        max_points_in_quadrant = num_points_per_sphere  # No restringir por defecto

    all_points = []
    radius_values = np.linspace(distance_start, distance_end, num_spheres)

    for radius in radius_values:
        valid_layer = False
        intentos = 0
        max_intentos = 100
        while not valid_layer and intentos < max_intentos:
            points = np.random.normal(size=(num_points_per_sphere, 3))
            norms = np.linalg.norm(points, axis=1)
            normalized_points = points / norms[:, None] * radius
            translated_points = normalized_points + center

            # Filtrar puntos en el cuadrante positivo
            filtered_points = translated_points[
                (translated_points[:, 0] > 0) &
                (translated_points[:, 1] > 0) &
                (translated_points[:, 2] > 0)]

            intentos += 1
            # Verificar si cumple con el rango de puntos en el cuadrante
            if min_points_in_quadrant <= filtered_points.shape[0] <= max_points_in_quadrant:
                valid_layer = True
                all_points.append(filtered_points)
            if intentos >= max_intentos:
                print(f"Advertencia: Se alcanzó el máximo de intentos ({max_intentos}) sin cumplir los criterios.")
                break

    # Combinar todos los puntos de todas las esferas en un solo array y redondear
    all_points = np.vstack(all_points) if all_points else np.array([])
    all_points = np.around(all_points, decimals=8)

    return all_points


def generate_3d_points_on_x(x_min, x_max, num_points):
    # Generar puntos uniformemente espaciados en el eje x
    x_points = np.linspace(x_min, x_max, num_points)

    # Crear arrays de ceros para y y z con el mismo número de puntos que x
    y_points = np.zeros(num_points)
    z_points = np.zeros(num_points)

    # Combinar x, y, y z en un único array 3D
    points_3d = np.vstack((x_points, y_points, z_points)).T  #Transponer para ajustar la forma

    return points_3d


def generate_rotation_combinations(nx, ny, nz):
    '''This function is to make the angle partition using the angle jump'''

    rotations = []
    for x in range(0, 360, nx):
        for y in range(0, 360, ny):
            for z in range(0, 360, nz):
                rotations.append((x, y, z))
    return rotations


def calcular_centroide_geometrico(puntos):
    """
    Calcula el centroide geométrico de un conjunto de puntos en 3D.
    """
    centroide = np.mean(puntos, axis=0)
    return centroide

def configurations_list(base_molecule, random_points, rotations=None, include_rotations=True):
    if rotations is None:
        rotations = []  # Establecer rotations a una lista vacía si no se proporciona
    
    centre_mol = calcular_centroide_geometrico(base_molecule)
    configurations = []
    for point in random_points:
        translation_vector = point - centre_mol
        translated_molecule = base_molecule + translation_vector
        if include_rotations and rotations:  # Verificar también que rotations no esté vacío
            for rotation in rotations:
                rot = R.from_euler('xyz', rotation, degrees=True)
                rotated_molecule = rot.apply(translated_molecule)
                configurations.append(np.array(rotated_molecule))
        else:
            configurations.append(np.array(translated_molecule))
    return configurations



ruta_archivo_xyz = r"C:\Users\migue\Downloads\CHN_mol_multipoles_miguel.xyz"
matriz_coordenadas = xyz_a_matriz_numpy(ruta_archivo_xyz)
print(matriz_coordenadas)


CHN_xaxis_rot_config = configurations_list(matriz_coordenadas,
                                           generate_3d_points_on_x(6, 9, 8),
                                           generate_rotation_combinations(180, 180, 180))

CHN_xaxis_rot_config = np.around(CHN_xaxis_rot_config, decimals=8)
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
name_vec = f'CHN_xaxis_rot_config_{timestamp}.npy'
np.save(name_vec, CHN_xaxis_rot_config)

