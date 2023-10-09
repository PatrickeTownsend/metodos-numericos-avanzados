import numpy as np
A = np.array([[0,2,3],[4,0,6],[7,8,0]])
B = [1,2,3]
pos = np.argmax(A[:,0])
pos2 = np.where(A == 0)
print(pos2[0])