# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 18:02:01 2024

@author: migue
"""

import psi4
import numpy as np
from datetime import datetime
import os

# Configura la precisión de impresión para NumPy
np.set_printoptions(precision=5, linewidth=200, threshold=2000, suppress=True)

# Obtén el nombre del script sin la extensión .py
script_name = os.path.splitext(os.path.basename(__file__))[0]


def set_up_psi4(timestamp):
    psi4.set_memory('200 GB')
    # Incorpora el nombre del script en el nombre del archivo de salida
    output_filename = f'{script_name}_{timestamp}.dat'
    psi4.core.set_output_file(output_filename, True)
    psi4.set_options({
        'basis': 'aug-cc-pvtz',
        'freeze_core': True,
        'e_convergence': 1e-10,
        'cc_num_threads': 12,
        'd_convergence': 1e-10,
        'ints_tolerance': 1e-18
    })

def leer_xyz_con_carga_multiplicidad_tipos_y_coordenadas(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    carga_multiplicidad = lines[0].strip()
    atomos_y_coordenadas = "\n".join(lines[3:])  # Usa directamente las líneas con átomos y coordenadas
    return carga_multiplicidad, atomos_y_coordenadas


def generar_geometria_psi4(archivo_mol1, archivo_mol2, coor):
    carga_multiplicidad_mol1, atomos_y_coordenadas_mol1 = leer_xyz_con_carga_multiplicidad_tipos_y_coordenadas(
        archivo_mol1)
    carga_multiplicidad_mol2, atomos_y_coordenadas_mol2 = leer_xyz_con_carga_multiplicidad_tipos_y_coordenadas(
        archivo_mol2)

    geometrias = []
    for config in coor:
        # Filtra las líneas vacías antes de enumerar para la molécula 2
        lineas_no_vacias_mol2 = [line for line in atomos_y_coordenadas_mol2.splitlines() if line.strip()]

        coordenadas_mol2 = "\n".join(
            f"{line.split()[0]:<2} {config[i, 0]:>12.8f} {config[i, 1]:>12.8f} {config[i, 2]:>12.8f}"
            for i, line in enumerate(lineas_no_vacias_mol2)
        )

        # Elimina líneas vacías entre las coordenadas de la primera molécula y ajusta el formato
        atomos_y_coordenadas_mol1_sin_espacios = "\n".join(
            f"{line.split()[0]:<2} {float(line.split()[1]):>12.8f} {float(line.split()[2]):>12.8f} {float(line.split()[3]):>12.8f}"
            for line in atomos_y_coordenadas_mol1.splitlines() if line.strip())

        geometria = f"""{carga_multiplicidad_mol1}
{atomos_y_coordenadas_mol1_sin_espacios}
--
{carga_multiplicidad_mol2}
{coordenadas_mol2}

units angstrom
symmetry c1
"""
        geometrias.append(geometria.strip())
    return geometrias


def calculate_energies(geometrias):
    esapt, eind, edisp, eelst, eexch = np.zeros((5, len(geometrias)))
    for i, geometria in enumerate(geometrias):
        psi4.geometry(geometria)
        psi4.energy('sapt2')
        esapt[i], eind[i], edisp[i], eelst[i], eexch[i] = [
            psi4.variable(name) * 627.509 for name in
            ('SAPT TOTAL ENERGY', 'SAPT IND ENERGY', 'SAPT DISP ENERGY', 'SAPT ELST ENERGY', 'SAPT EXCH ENERGY')
        ]
        psi4.core.clean()
    return esapt, eind, edisp, eelst, eexch


def save_results_npz(esapt, eind, edisp, eelst, eexch, filename):
    np.savez(filename, esapt=esapt, eind=eind, edisp=edisp, eelst=eelst, eexch=eexch)

# Ejecución principal
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
set_up_psi4(timestamp)

# Aquí necesitas cargar tus coordenadas y definir las rutas a tus archivos XYZ
coor = np.load('/home/miguel/QM_Project/pf6_quin_energies_sapt2/Pf6_quin_spherical_config_20240226_110016.npy')
archivo_mol1 = '/home/miguel/QM_Project/pf6_quin_energies_sapt2/PF6_centered_withcharge.xyz'
archivo_mol2 = '/home/miguel/QM_Project/pf6_quin_energies_sapt2/CHN_multipole_charged.xyz'

geometrias = generar_geometria_psi4(archivo_mol1, archivo_mol2, coor)
energies = calculate_energies(geometrias)

filename = f'{script_name}_resultados_{timestamp}.npz'
save_results_npz(*energies, filename=filename)
# remember that different input files should have different configurations file
