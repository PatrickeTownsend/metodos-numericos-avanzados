import numpy as np
import functions as func
import solvers as solv
import matplotlib.pyplot as plt
N = 10000
dt = 0.001
x = np.zeros((N,4))
y = np.zeros((N,4))
z = np.zeros((N,4))
att = [func.Lorenz_atractor,func.Rossler_attractor,func.Thomas_attractor,func.Rossler_attractor]
for j in range(4):
  t = 0.0
  for i in range(N-1):
    x[:,j],y[:,j],z[:,j] = solv.Runge_Kutta4(t,x[:,j],y[:,j],z[:,j],att[j],dt,i)
    t += dt

plt.figure()
ax = plt.axes(projection ='3d')
ax.plot3D(x[:,0],y[:,0],z[:,0])
plt.show()


