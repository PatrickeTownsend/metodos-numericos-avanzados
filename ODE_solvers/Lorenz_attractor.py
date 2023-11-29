import numpy as np
import solvers as solv
import matplotlib.pyplot as plt
def f1(t,x,y,z):
    dx = 10*(-x+y)
    dy = -x*z + 28*x - y
    dz = x*y -(8/3)*z
    return dx,dy,dz

def Lorenz(x0,y0,z0,dt,N):
    x = np.zeros((N,3))
    y = np.zeros((N,3))
    z = np.zeros((N,3))
    t = 0.0
    for i in range(N-1):
        x[:,0],y[:,0],z[:,0] = solv.Runge_Kutta4(t,x[:,0],y[:,0],z[:,0],f1,dt,i)
        x[:,1],y[:,1],z[:,1] = solv.Euler_exp(t,x[:,0],y[:,0],z[:,0],f1,dt,i)
    
    fig, (ax1,ax2) = plt.subplots(1,2)
    fig.suptitle("Lorenz attractor")
    axi = plt.axes(projection="3d")
    ax1.plot3D(x[:,0],y[:,0],z[:,0])
    ax1.title("Runge-Kutta 4")
    ax2.plot3D(x[:,1],y[:,1],z[:,1])
    ax2.title("Euler explícito")
    plt.show()
    




