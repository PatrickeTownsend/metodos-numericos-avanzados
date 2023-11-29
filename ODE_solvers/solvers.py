import functions as func

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

