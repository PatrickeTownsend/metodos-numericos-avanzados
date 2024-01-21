import math
import solvers as solv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
def Lorenz_attractor(t,x,y,z):
    dx = 10*(-x+y)
    dy = -x*z + 28*x - y
    dz = x*y -(8/3)*z
    return dx,dy,dz

def Rossler_attractor(t,x,y,z):
    dx = -(y+z)
    dy = x + 0.2*y
    dz = 0.2 + z*(x-5.7)
    return dx,dy,dz

def Thomas_attractor(t,x,y,z):
    dx = math.sin(y) - 0.208186*x
    dy = math.sin(z) - 0.208186*y
    dz = math.sin(x) - 0.208186*z
    return dx,dy,dz

def Three_Scroll_attractor(t,x,y,z):
    dx = 32.48*(y-x) + 0.13*x*z
    dy = 45.84*x - x*z + 14.7*y
    dz = 1.18*z + x*y - 0.57*x**2
    return dx,dy,dz

def Lorenz_J(x,y,z,dt):
    J = np.array([[-1-10*dt, 10*dt, 0],
                  [-dt*z+28*dt, -1-dt, -dt*x],
                  [dt*y, dt*x, -1-(8/3)*dt]],dtype=float)
    return J
def Rossler_J(x,y,z,dt):
    J = np.array([[-1, -dt, dt],
                  [dt, -1+0.2*dt, 0],
                  [dt*z, 0, -1+dt*x-5.7*dt]])
    return J
def Thomas_J(x,y,z,dt):
    J = np.array([[-1-dt*0.208186, dt*np.cos(y), 0],
                  [0, -1-dt*0.208186, dt*np.cos(z)],
                  [dt*np.sin(x), 0, -1-dt*0.208186]])
    return J
def ThreeScoll_J(x,y,z,dt):
    J = np.array([[-1-dt*32.48+0.13*dt*z, dt*32.48, 0.13*dt*x],
                  [45.84*dt-dt*z, -1+14.7*dt, -dt*x],
                  [dt*y-2*0.57*x*dt, dt*x, -1+1.18*dt]])
    return J

def Lorenz(x0,y0,z0,dt,N):
    x = np.zeros((N,3))
    y = np.zeros((N,3))
    z = np.zeros((N,3))

    x[0,:] = x0
    y[0,:] = y0
    z[0,:] = z0
    t = 0.0
    for i in range(N-1):
        x[:,0],y[:,0],z[:,0] = solv.Runge_Kutta4(t,x[:,0],y[:,0],z[:,0],Lorenz_attractor,dt,i)
        x[:,1],y[:,1],z[:,1] = solv.Euler_exp(t,x[:,1],y[:,1],z[:,1],Lorenz_attractor,dt,i)
        x[:,2],y[:,2],z[:,2] = solv.Euler_imp(t,x[:,2],y[:,2],z[:,2],Lorenz_attractor,dt,i,Lorenz_J)
    
    fig, ax1 = plt.subplots(1, 3, figsize=(10, 5), subplot_kw={'projection': '3d'})


    ax1[0].plot(x[:,0], y[:,0], z[:,0])
    ax1[0].set_title('Runge-Kutta 4')
    ax1[0].set_xlabel("x")
    ax1[0].set_ylabel("y")
    ax1[0].set_zlabel("z")


    ax1[1].plot(x[:,1], y[:,1], z[:,1])
    ax1[1].set_title('Euler explicito')
    ax1[1].set_xlabel("x")
    ax1[1].set_ylabel("y")
    ax1[1].set_zlabel("z")

    ax1[2].plot(x[:,2], y[:,2], z[:,2])
    ax1[2].set_title('Euler implicito')
    ax1[2].set_xlabel("x")
    ax1[2].set_ylabel("y")
    ax1[2].set_zlabel("z")
    fig.suptitle("Atractor de Lorenz")
    #plt.savefig("ODE_solvers/plots/lorenz.png")
    plt.show()
    
def Thomas(x0,y0,z0,dt,N):
    x = np.zeros((N,3))
    y = np.zeros((N,3))
    z = np.zeros((N,3))

    x[0,:] = x0
    y[0,:] = y0
    z[0,:] = z0
    t = 0.0
    for i in range(N-1):
        x[:,0],y[:,0],z[:,0] = solv.Runge_Kutta4(t,x[:,0],y[:,0],z[:,0],Thomas_attractor,dt,i)
        x[:,1],y[:,1],z[:,1] = solv.Euler_exp(t,x[:,1],y[:,1],z[:,1],Thomas_attractor,dt,i)
        x[:,2],y[:,2],z[:,2] = solv.Euler_imp(t,x[:,2],y[:,2],z[:,2],Thomas_attractor,dt,i,Thomas_J)
    
    fig, ax1 = plt.subplots(1, 3, figsize=(10, 5), subplot_kw={'projection': '3d'})


    ax1[0].plot(x[:,0], y[:,0], z[:,0])
    ax1[0].set_title('Runge-Kutta 4')
    ax1[0].set_xlabel("x")
    ax1[0].set_ylabel("y")
    ax1[0].set_zlabel("z")


    ax1[1].plot(x[:,1], y[:,1], z[:,1])
    ax1[1].set_title('Euler explicito')
    ax1[1].set_xlabel("x")
    ax1[1].set_ylabel("y")
    ax1[1].set_zlabel("z")

    ax1[2].plot(x[:,2], y[:,2], z[:,2])
    ax1[2].set_title('Euler implicito')
    ax1[2].set_xlabel("x")
    ax1[2].set_ylabel("y")
    ax1[2].set_zlabel("z")
    fig.suptitle("Atractor de Thomas")
    #plt.savefig("ODE_solvers/plots/thomas.png")
    plt.show()
    
def Rossler(x0,y0,z0,dt,N):
    x = np.zeros((N,3))
    y = np.zeros((N,3))
    z = np.zeros((N,3))

    x[0,:] = x0
    y[0,:] = y0
    z[0,:] = z0
    t = 0.0
    for i in range(N-1):
        x[:,0],y[:,0],z[:,0] = solv.Runge_Kutta4(t,x[:,0],y[:,0],z[:,0],Rossler_attractor,dt,i)
        x[:,1],y[:,1],z[:,1] = solv.Euler_exp(t,x[:,1],y[:,1],z[:,1],Rossler_attractor,dt,i)
        x[:,2],y[:,2],z[:,2] = solv.Euler_imp(t,x[:,2],y[:,2],z[:,2],Rossler_attractor,dt,i,Rossler_J)
    
    fig, ax1 = plt.subplots(1, 3, figsize=(10, 5), subplot_kw={'projection': '3d'})


    ax1[0].plot(x[:,0], y[:,0], z[:,0])
    ax1[0].set_title('Runge-Kutta 4')
    ax1[0].set_xlabel("x")
    ax1[0].set_ylabel("y")
    ax1[0].set_zlabel("z")


    ax1[1].plot(x[:,1], y[:,1], z[:,1])
    ax1[1].set_title('Euler explicito')
    ax1[1].set_xlabel("x")
    ax1[1].set_ylabel("y")
    ax1[1].set_zlabel("z")

    ax1[2].plot(x[:,2], y[:,2], z[:,2])
    ax1[2].set_title('Euler implicito')
    ax1[2].set_xlabel("x")
    ax1[2].set_ylabel("y")
    ax1[2].set_zlabel("z")
    fig.suptitle("Atractor de Rossler")
    #plt.savefig("ODE_solvers/plots/rossler.png")
    plt.show()
    

def Three_Scroll(x0,y0,z0,dt,N):
    x = np.zeros((N,3))
    y = np.zeros((N,3))
    z = np.zeros((N,3))

    x[0,:] = x0
    y[0,:] = y0
    z[0,:] = z0
    t = 0.0
    for i in range(N-1):
        x[:,0],y[:,0],z[:,0] = solv.Runge_Kutta4(t,x[:,0],y[:,0],z[:,0],Three_Scroll_attractor,dt,i)
        x[:,1],y[:,1],z[:,1] = solv.Euler_exp(t,x[:,1],y[:,1],z[:,1],Three_Scroll_attractor,dt,i)
        x[:,2],y[:,2],z[:,2] = solv.Euler_imp(t,x[:,2],y[:,2],z[:,2],Three_Scroll_attractor,dt,i,ThreeScoll_J)
    
    fig, ax1 = plt.subplots(1, 3, figsize=(10, 5), subplot_kw={'projection': '3d'})


    ax1[0].plot(x[:,0], y[:,0], z[:,0])
    ax1[0].set_title('Runge-Kutta 4')
    ax1[0].set_xlabel("x")
    ax1[0].set_ylabel("y")
    ax1[0].set_zlabel("z")


    ax1[1].plot(x[:,1], y[:,1], z[:,1])
    ax1[1].set_title('Euler explicito')
    ax1[1].set_xlabel("x")
    ax1[1].set_ylabel("y")
    ax1[1].set_zlabel("z")

    ax1[2].plot(x[:,2], y[:,2], z[:,2])
    ax1[2].set_title('Euler implicito')
    ax1[2].set_xlabel("x")
    ax1[2].set_ylabel("y")
    ax1[2].set_zlabel("z")
    fig.suptitle("Atractor de Three-Scroll")
    #plt.savefig("ODE_solvers/plots/three_scroll.png")
    plt.show()
    
