from pyEnGNet.Notebooks.Code.pyEnGNet import PyEnGNet
import numpy as np
import time as t

fils = 500
numcores = 8


print("EN PROGRESO: TESTS DE TIEMPOS PARA MATRIZ 500 X 500")
with open("../Data/Performance measures/t_num_cores_500x500.csv", "a") as f:
    f.write("numcores,tiempo(s)\n")
    for ncores in [1, 2, 4, 6, 8]:
        print("Calculando con", ncores, "cores")
        maindata = np.random.random_sample((fils, fils))
        peg = PyEnGNet(maindata, cores=ncores)
        t0 = t.time()
        peg.engnet_1_0()
        tf = t.time() - t0
        print("Tiempo total:", tf)
        f.write(str(ncores)+","+str(tf)+"\n")
print("FINALIZADO: TESTS DE TIEMPOS PARA MATRIZ 500 X 500")
'''
    for filas in [10, 50, 100, 500, 1000, 2000, 5000]:
        maindata = np.random.random_sample((filas, cols))
        peg = PyEnGNet(maindata, cores=1)
        peg.engnet_1_0()
        
    '''