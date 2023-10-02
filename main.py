import numpy as np
def suma(x,y):
    z = x+y
    return z
def LUdcmp(A):
    N = np.shape(a)[0]
    for j in range(1,N+1):
        for i in range(j+1,N+1):
            l = U[i+1,j+1] / U[j+1,j+1]
            
    return L,U





if __name__ == "__main__":
    
    a = np.array([[1,2,3],[4,6,5],[7,8,7]])
    b = np.array([1,2,3])
    L, U = LUdcmp

    c = suma(a,b)
    print("c",c)
    