import functions as func
import numpy as np
def Runge_Kutta4(t:float,xi,yi,zi,f,dt:float,i):
    x = xi[i]
    y = yi[i]
    z = zi[i]

    k1 = f(t,x,y,z)
    k2 = f(t + dt/2.0, x + k1[0]*dt/2.0, y + k1[1]*dt/2.0, z + k1[2]*dt/2.0)
    k3 = f(t + dt/2.0, x + k2[0]*dt/2.0, y + k2[1]*dt/2.0, z + k2[2]*dt/2.0)
    k4 = f(t + dt, x + k3[0]*dt, y + k3[1]*dt, z + k3[2]*dt)

    xi[i+1] = x + (dt/6)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
    yi[i+1] = y + (dt/6)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
    zi[i+1] = z + (dt/6)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
    return xi,yi,zi

def Euler_exp(t,x,y,z,f,dt,i):
    x[i+1] = x[i] + dt*f(t,x[i],y[i],z[i])[0]
    y[i+1] = y[i] + dt*f(t,x[i],y[i],z[i])[1]
    z[i+1] = z[i] + dt*f(t,x[i],y[i],z[i])[2]

    return x,y,z

def Newthon_Raphson(t,x0,y0,z0,dt,f,Jac,max_tol,max_iter):
    x = x0
    y = y0
    z = z0
    for i in range(max_iter):
        fx = x0 + dt*f(t,x,y,z)[0] - x
        fy = y0 + dt*f(t,x,y,z)[1] - y
        fz = z0 + dt*f(t,x,y,z)[2] - z

        f_new = np.array([[fx],[fy],[fz]])
        J = Jac(x,y,z,dt)
        delta = np.linalg.solve(J,f_new)
        x = x - delta[0,0]
        y = y - delta[1,0]
        z = z - delta[2,0]

        if np.linalg.norm(delta) < max_tol:
            return x,y,z

def Euler_imp(t,x,y,z,f,dt,i,Jac):
    x[i+1],y[i+1],z[i+1] = Newthon_Raphson(t,x[i],y[i],z[i],dt,f,Jac,1e-8,1000)
    return x,y,z
    
         