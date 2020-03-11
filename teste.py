import numpy as np

rigidez = np.array([1.59, -0.40,-0.54])

forcas = np.array([0, 150, -100])

print(rigidez - forcas)