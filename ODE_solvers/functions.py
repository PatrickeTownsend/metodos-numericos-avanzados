import math
def Lorenz_atractor(t,x,y,z):
    dx = 10*(-x+y)
    dy = -x*z + 28*x - y
    dz = x*y -(8/3)*z
    print(dy)
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

def Three_Scroll(t,x,y,z):
    dx = 32.48*(y-x) + 0.13*x*z
    dy = 45.84*x - x*z + 14.7*y
    dz = 1.18*z + x*y - 0.57*x**2
    return dx,dy,dz

def init_lorenz(x,y,z):
    x[0] = 1.1
    y[0] = 2.0
    z[0] = 7.0
    return x,y,z
