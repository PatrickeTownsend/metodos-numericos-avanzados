import numpy as np

def LUdcmp(A):
    U = A
    L = np.eye(np.shape(A)[0])
    N = np.shape(A)[0]
    for j in range(1,N+1):
        for i in range(j+1,N+1):
            L[i-1,j-1] = U[i-1,j-1] / U[j-1,j-1]
            U[i-1,:] =   U[i-1,:] - L[i-1,j-1]*U[j-1,:]
    return L,U


def BackSub(L, U, b):
    N = np.shape(L)[0]
    y = np.zeros((N,1))
    x = np.zeros((N,1))
    # Forward substitution
    y[0] = b[0]
    for i in range(1,N+1):
        y[i-1] = b[i-1] 
        for j in range(1,i):
            y[i-1] = y[i-1] - y[j-1]*L[i-1,j-1]
    # backward substitution
    x[N-1] = y[N-1]
    for i in range(N,0,-1):
        x[i-1] = y[i-1]
        for j in range(i+1,N+1):
            x[i-1] = x[i-1] - x[j-1]*U[i-1,j-1]
        x[i-1] = x[i-1] / U[i-1,i-1]
    return x

def Singular(A):
    (row,col) = np.shape(A)
    for i in range(0,row):
        for j in range(0,col):
            if i==j:
                if A[i,j]==0:
                    singular = True
    return singular

if __name__ == "__main__":
    
    A = np.array([[1,2,3],[4,6,5],[7,8,7]])
    b = np.array([1,1,2])
    
    if Singular(A):
        for j in range(0,np.shape(A)[1]):
            col = list(A[:,j])
            zero_pos = col.index(0)

    L, U = LUdcmp(A)
    x = BackSub(L,U,b)

    print("L",L)
    print("\n")
    print("U",U)
    print("\n")
    print("A:", np.dot(L,U))
    print("\n")
    print("x:",x)
    