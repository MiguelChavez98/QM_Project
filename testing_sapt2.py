import psi4 
import numpy as np 
import scipy 
np.set_printoptions(precision=5, linewidth=200, threshold=2000, suppress=True)


numpy_memory = 100

psi4.set_memory('100 GB')

memory 100 GB

psi4.core.set_output_file('output_chn_sapt2_v2.dat', True)

psi4.set_options({'basis': 'aug-cc-pvtz',
                  'freeze_core': True,
                  'e_convergence': 1e-10,
                  'cc_num_threads': 12,
                  'd_convergence': 1e-10,
                  'ints_tolerance': 1e-18})

distances = [5.7]
esapt = np.zeros((len(distances)))
eind = np.zeros((len(distances)))
edisp = np.zeros((len(distances)))
eelst = np.zeros((len(distances)))
eexch = np.zeros((len(distances)))

for i in range(len(distances)):
  dimer = psi4.geometry("""
  1 1
C    0.04823180   0.00000000   1.29854000
C    0.77157200  -1.24836000   0.73279000
C    0.77157200   1.24836000   0.73279000
C   -1.43595000   0.00000000  -0.70655000
C    0.68327200  -1.24567000  -0.76739000
C    0.68327200   1.24567000  -0.76739000
C   -1.40934000   0.00000000   0.78001000
H   -0.08271820   0.00000100  -2.22938000
H    0.07010180  -0.00000100   2.29838000
H    1.72063000  -1.24122000   1.01304000
H    1.72063000   1.24122000   1.01304000
H    0.34934200  -2.06980000   1.09097000
H    0.34934200   2.06980000   1.09097000
H   -1.94519000   0.80829000  -1.08795000
H   -1.94519000  -0.80829000  -1.08795000
H    0.20046200  -2.05221000  -1.07690000
H    0.20046200   2.05221000  -1.07690000
H    1.59420000  -1.26188000  -1.15623000
H    1.59420000   1.26188000  -1.15623000
H   -1.94519000  -0.88912000   1.14598000
H   -1.94519000   0.88912000   1.14598000
N   -0.04853820   0.00000000  -1.22961000
--
1 1
C    6.04823341   0.00000000   1.29853655
C    6.77157441  -1.24835900   0.73279055
C    6.77157441   1.24835900   0.73279055
C    4.56405441   0.00000000  -0.70654845
C    6.68327341  -1.24566500  -0.76739145
C    6.68327341   1.24566500  -0.76739145
C    4.59066541   0.00000000   0.78001155
H    5.91728241   0.00000900  -2.22937745
H    6.07010241  -0.00000900   2.29838155
H    7.72062641  -1.24121900   1.01303855
H    7.72062641   1.24121900   1.01303855
H    6.34933941  -2.06979700   1.09097155
H    6.34933941   2.06979700   1.09097155
H    4.05481241   0.80829000  -1.08795045
H    4.05481241  -0.80829000  -1.08795045
H    6.20046141  -2.05221200  -1.07689845
H    6.20046141   2.05221200  -1.07689845
H    7.59419841  -1.26188400  -1.15623045
H    7.59419841   1.26188400  -1.15623045
H    4.05481241  -0.88911900   1.14597555
H    4.05481241   0.88911900   1.14597555
N    5.95146541   0.00000000  -1.22961445

units angstrom
symmetry c1
""")

  psi4.energy('sapt2', molecule=dimer)
  esapt[i] = psi4.variable('SAPT TOTAL ENERGY')*627.509
  eind[i] = psi4.variable('SAPT IND ENERGY')*627.509
  edisp[i] = psi4.variable('SAPT DISP ENERGY')*627.509
  eelst[i] = psi4.variable('SAPT ELST ENERGY')*627.509
  eexch[i] = psi4.variable('SAPT EXCH ENERGY')*627.509
  psi4.core.clean()


# Convertir el vector en una cadena de texto
vector_str = str(esapt)
# Abrir un archivo en modo de escritura (esto creará el archivo si no existe)
with open('sapt2_quin.txt', 'w') as archivo:
    # Escribir la cadena en el archivo
    archivo.write(vector_str)


vector_d_str = str(distances)
with open('sapt2_quin.txt', 'a') as archivo:
    # Es posible que quieras añadir un salto de línea antes del nuevo contenido
    archivo.write('\n')  # Esto añade una nueva línea
    # Escribir la nueva cadena en el archivo
    archivo.write(vector_d_str)
print(esapt)


vector_d_str = str(edisp)
with open('sapt2_quin.txt', 'a') as archivo:
    # Es posible que quieras añadir un salto de línea antes del nuevo contenido
    archivo.write('\n')  # Esto añade una nueva línea
    # Escribir la nueva cadena en el archivo
    archivo.write(vector_d_str)
print(esapt)


vector_d_str = str(eind)
with open('sapt2_quin.txt', 'a') as archivo:
    # Es posible que quieras añadir un salto de línea antes del nuevo contenido
    archivo.write('\n')  # Esto añade una nueva línea
    # Escribir la nueva cadena en el archivo
    archivo.write(vector_d_str)
print(esapt)


vector_d_str = str(eelst)
with open('sapt2_quin.txt', 'a') as archivo:
    # Es posible que quieras añadir un salto de línea antes del nuevo contenido
    archivo.write('\n')  # Esto añade una nueva línea
    # Escribir la nueva cadena en el archivo
    archivo.write(vector_d_str)
print(esapt)


vector_d_str = str(eexch)
with open('sapt2_quin.txt', 'a') as archivo:
    # Es posible que quieras añadir un salto de línea antes del nuevo contenido
    archivo.write('\n')  # Esto añade una nueva línea
    # Escribir la nueva cadena en el archivo
    archivo.write(vector_d_str)
print(esapt)
