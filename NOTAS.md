Notas para el proyecto actualizaciones etc etc. 

# May 4 2024

Hoy toca hacer una optimizacion de los parametros que habia obtenido. So one of the first things to do is check all the data that is availible in the github repository and then judge what points are we going to use and the optimization methods. 

We have the information of 3 configurations of PF6 and the energy of the dimers that arise from the combination of that configurations and 19 rotations of the molecule quin. 

so we are going to take all these 19X3 values and try to fit them first to a new set of values A B and C for the molecules. 

We also have to check the order of the molecules and the values that are being used. 

The XYZ files are in angstroms, OpenMM uses the values of the molecules in angstroms because the PDB file is already consider in angstroms.

I need to check if all the values that are being used in the files are angstroms or maybe they are nm o some kind of different measure system like bohr ? because that is what I used in the psi4 library, in the code I know that I specified that the distances were angstroms but need to check the original files to make sure of that. 


I have the original files in the github project and also I have the first first moelcule model that alston send me with the geometry corrected. 


After that I can assure that the molecules are correct and that I need to check the variables and units of the constants. 

-- All seems to be in order with the archives(checked in mercury) the archives that I used and the ones that anthony sendme have the same lengths, I have some archives with different lengths but I dont know why they have different lengths, but that ones were not used, maybe they are the right ones? idk but at least im using the ones that anthony send me. Actually can they ( the bad ones) be realted with the energy minimization that alston or anthony did in some time ago? who knows.

# May 6 2024 

Now that I have have all the coordinates right I will consider the units of the values that I have for the constants. 

**values for the dispersion constants from camcasp (Alston)**

Para las constantes de alston los valores son correctos sin embargo las transformaciones eran solo numericas me parece que no me tenia que meter con las unidades mas alla de los valores que se ponen en la formula. 

from the fitting for FP6 

Alston's values 

Dispercion P = 30.61343 hartree bohr**6  # factor necesario es 13.77929292 y la k 
Dispersion F = 11.50171 hartree bohr**6  # factor necesario es 13.77929292 y la k

Ahora los obtenidos por medio de la optimizacion: 

BF = 305.98575047049286 # que serian AB en orden de la formula #A  
CF = 3.2512195680784135 # que serian AB en el orden de la formula #B

BP = 48058.29804395409 # que serian AB en orden de la formula #A
CP = 1.72034397464677  # que serian AB en el orden de la formula #B

la constante para escalar la dispersion es la k deeee 

K = -5.203428163043108 # la k esta negativa por que deacuerdo al fitting como se suma todo al final pero se supone que se reste la dispersion pues es negativa. 

# May 7 2024

Ahora los valores encontrados son con constantes normales pero pero pero, Como todas las formulas usan la constante multiplicando por bohr las distancias y transformando a kcal lo importante es adecuar las constantes con estos terminos ya incluidos en todo lo que se esta haciendo. 

Los BF y BP y CF y CP de las optimizaciones son constantes que estan siempre afectadas por transformaciones antes de salir y por eso mismo creo conveniente agregar esas transformaciones en las constantes si es que las vamos a usar con openmm

el factor para la C por ejemplo es solamente aplicar angstrom a bohr y se vuelve en el factor

Para el caso de la B es solamente aplicar La conversion final pues es una multiplicacion cualquiera de todos los terminos por lo que puedes ponerlo ahi metido la conversion de hartree a kcal 

Ahora vamos a ver, para no convertir a kcalorias sabemos que una kilo caloria es 4.184 kJ/mol

Entonces vamos a utilizar los factores para llevar las constantes a eso 

Para el caso de las [B] Es necesario tener que multiplicar al final por 627.509 y luego por 4.184

Dando un factor de [2625.497656]

Para el caso de la constante [C] Lo que tenemos es que convertir a bohr primero asi que lo asignamos desde el inicio: 

Factor de Angstrom to bohr es: [1.8897261246257702]

Ahora para el caso de la dispersion lo que tenemos es una valor de k de 5.203428163043108 y multiplicando eso por el valor de el factor de 13.77929292 queda un factor para esas cosntantes de: 
[71.699560] Esta constante es de kcalorias por que en realidad lo que tenemos son kilocalorias como factor y la parte de bohr a la 6 ya se hizo para poderlo incluir en la constante de dispersion y si le agregamos que ponemos tambien el valor de conversion pues ya casi esta todo en esa constante pero como es kcal lo pasamos a kj con el 4. 

Esto esta mal, para el caso de la dispersion lo que tenemos es una constante que si ponemos todo en unidades atomicas nos deberia dar hartree por que es lo que las unidades originales tienen, asi queajustamos los valores para obtener unidades atomicas y despues podemos poner la conversion a hartree dentro de las misma variable o constante, asi que ademas de eso hay que ver como se eliminan y multiplican los datos por que si lo que buscamos es tener hartree multiplicar por 13 que es el factor de angstrom a bohr hace que el resultado ya no este en hartree, mas bien lo tenemos que incluir en las constantes desde ahorita. #may 10 2024

seriaaaa, mentira si lo hice bien pero de menso no lo escribi, por que tengo el 1.8897 a la 6 abajo, bien, pero tengo arriba el de convertir a kcal asi que ahi realmente no esta aca, mas bien el que falta ahi es el por 4 por que son kilocalorias lo que estoy convirtiendo el hartree y la constante a si que vamos a checar eso # may 10 2024

Ahora ya tenemos los 3 factores para las 3 constantes. Vamos a escribirlas y poner los factores: 


Ap = 30.61343 * 71.699560 * 4.184 = 9183.752225203909 ESte seria C Este deberia ser correcto--
ESTOS PODRIAN IR ENTRE 2 osea divididos. 
Af = 11.50171 * 71.699560 * 4.184 = 3450.4090134999587

Bp = 48058.29804395409 * 2625.497656 = 126176948.86575085 # :CCC Este seria A
Bf = 305.98575047049286 * 2625.497656 = 803364.8706296799 # :/

Cp = 1.72034397464677 * 1.8897261246257702 = 3.2509789522325345      Este seria B 
Cf = 3.2512195680784135 * 1.8897261246257702 = 6.143914554692291

# May 12

Archivo de penalizacion llamado penalizacion.py en SPIDER para aplicar a la parte de la optimizacion. 

Basicamente vamos a crear un vector de penalizaciones para tener un vector que sirva para castigar de forma individual que tanto se aleja un parametro de lo esperado. 



# May 14 

Orientation time correlation vector dot poduct in time zero. Position

Lower temperatures

box around. 

super cell. 

automatization. 

try models like argon in crystals like fcc. 

gulp models 

----


send you the model as it is --|

wringtng 

input file pdb file. 

