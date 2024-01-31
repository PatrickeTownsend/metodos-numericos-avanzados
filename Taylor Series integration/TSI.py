import numpy as np
import matplotlib.pyplot as plt

# Constantes
mu_cb = 3.986004418 * 1e5

# Parámetros de la integración

timeinit = 0
timestop = 3.821*24*3600
timestep = 10
t = np.linspace(timeinit, timestop, int((timestop-timeinit)/timestep)+1)
steps = round((timestop - timeinit) / timestep)
k = 20
T = 2e-4
m0 = 100
Isp = 3000
g0 = 9.80665
mdot = -T/(Isp*g0) 


# Condiciones Iniciales

x = np.zeros((11,steps+1))
x[0,0] = 6678
x[1,0] = 0
x[2,0] = 0
x[3,0] = 0
x[4,0] = np.sqrt(mu_cb / (np.sqrt(x[0,0]**2 + x[1,0]**2 + x[2,0]**2)))
x[5,0] = 0
x[6,0] = m0
x[7,0] = x[0,0]**2 + x[1,0]**2 + x[2,0]**2
x[8,0] = x[7,0]**(3/2)
x[9,0] = x[3,0]**2+x[4,0]**2+x[5,0]**2
x[10,0] = x[9,0]**(1/2) 



# Función u en el instante inicial
# def compute_u0(mu_cb, x):
#     u0 = np.array([x[3], x[4], x[5], -mu_cb*x[0]/x[8], -mu_cb*x[1]/x[8],
#                    -mu_cb*x[2]/x[8], 0, 2*x[0]*x[3]+2*x[1]*x[4]+2*x[2]*x[5],
#                    3/2*x[8]*(2*x[0]*x[3]+2*x[1]*x[4]+2*x[2]*x[5])/x[7]])
#     return u0
def compute_u0(mu_cb, x):
    vmod = ((x[3]**2+x[4]**2+x[5]**2)**(1/2))
    u0 = np.array([x[3], x[4], x[5], -mu_cb*x[0]/x[8]+T*x[3]/x[6]/vmod,
                                     -mu_cb*x[1]/x[8]+T*x[4]/x[6]/vmod,
                                     -mu_cb*x[2]/x[8]+T*x[5]/x[6]/vmod,
                                     mdot, 2*x[0]*x[3]+2*x[1]*x[4]+2*x[2]*x[5],
                                     3/2*x[8]*(2*x[0]*x[3]+2*x[1]*x[4]+2*x[2]*x[5])/x[7],
                                     2*x[3]*(-mu_cb*x[0]/x[8]+T*x[3]/x[6]/vmod)+2*x[4]*(-mu_cb*x[1]/x[8]+T*x[4]/x[6]/vmod)+2*x[5]*( -mu_cb*x[2]/x[8]+T*x[5]/x[6]/vmod),
                                     1/2*x[9]**(-1/2)
                                     ])
    return u0

# Calcula U en el paso inicial
u0 = compute_u0(mu_cb, x[:,0])
U = np.zeros((11,k+1))
U[:,0] = u0

# Calcula W (var. auxiliares) en el paso inicial

def compute_w0(U,x):
    w0 = np.array([
                    0,                  #W1
                    0,                  #W2
                    0,                  #W3
                    x[0]/x[8],          #W4
                    x[1]/x[8],          #W5
                    x[2]/x[8],          #W6
                    0,                  #W7
                    x[0]*x[3],          #W8
                    x[1]*x[4],          #W9
                    x[2]*x[5],          #W10
                    x[8]*U[7,0]/x[7],   #W11
                    x[3]/x[10],         #W12
                    x[4]/x[10],         #W13
                    x[5]/x[10],         #W14
                    x[3]*U[3,0],        #W15
                    x[4]*U[4,0],        #W16
                    x[5]*U[5,0],        #W17
                    x[9]**(-1/2),       #W18
                   ])
    return w0

W = np.zeros((18,k+1))
W[:,0] = compute_w0(U,x[:,0])



def recurrence(U,W,x,k_max,mu_cb):
    for k in range(1,k_max+1):
        term4 = 0
        term5 = 0
        term6 = 0
        term81 = 0
        term82 = 0
        term83 = 0
        term91 = 0
        term92 = 0
        term12 = 0
        term13 = 0
        term14 = 0
        term15 = 0
        term16 = 0
        term17 = 0
        for j in range(1,k+1):
            term4 += U[8,j-1]/j*W[3,k-j]
            term5 += U[8,j-1]/j*W[4,k-j]
            term6 += U[8,j-1]/j*W[5,k-j]
            term91 += U[7,j-1]/j*U[7,k-j]
            term92 += U[7,j-1]/j*W[8,k-j]
            term12 += U[10,j-1]/j*W[3,k-j]
            term13 += U[10,j-1]/j*W[4,k-j]
            term14 += U[10,j-1]/j*W[5,k-j]
        for j in range(1,k):
            term81 += U[0,j-1]/j*U[3,k-j-1]/(k-j)
            term82 += U[1,j-1]/j*U[4,k-j-1]/(k-j)
            term83 += U[2,j-1]/j*U[5,k-j-1]/(k-j)
        W[3,k] = 1/x[8]*(U[0,k-1]/k-term4)
        W[4,k] = 1/x[8]*(U[1,k-1]/k-term5)
        W[5,k] = 1/x[8]*(U[2,k-1]/k-term6)
        W[6,k] = 0
        W[7,k] = x[0]*U[3,k-1]/k+U[0,k-1]/k*x[3]+term81
        W[8,k] = x[1]*U[4,k-1]/k+U[1,k-1]/k*x[4]+term82
        W[9,k] = x[2]*U[5,k-1]/k+U[2,k-1]/k*x[5]+term83
        W[11,k] = 1/x[10]*(U[3,k-1]/k-term12)
        W[12,k] = 1/x[10]*(U[4,k-1]/k-term13)
        W[13,k] = 1/x[10]*(U[5,k-1]/k-term14)
        U[0,k] = U[3,k-1]/k
        U[1,k] = U[4,k-1]/k
        U[2,k] = U[5,k-1]/k
        U[3,k] = -mu_cb*W[3,k]+T/x[6]*W[11,k]
        for j in range(1,k):
            term15 += U[3,j-1]/j*U[3,k-j-1]
            term16 += U[4,j-1]/j*U[4,k-j-1]
            term17 += U[5,j-1]/j*U[5,k-j-1]
        W[14,k] = x[3]*U[3,k]+U[3,k-1]/k*U[3,k]+term15
        U[4,k] = -mu_cb*W[4,k]+T/x[6]*W[12,k]
        W[15,k] = x[4]*U[4,k]+U[4,k-1]/k*U[4,k]+term16
        U[5,k] = -mu_cb*W[5,k]+T/x[6]*W[13,k]
        W[16,k] = x[5]*U[5,k]+U[5,k-1]/k*U[5,k]+term17
        U[6,k] = 0
        W[10,k] = 1/x[7]*(x[8]*U[7,k]+term91-term92)
        U[7,k] = 2*W[7,k]+2*W[8,k]+2*W[9,k]
        U[8,k] = 3/2*W[10,k]
        U[9,k] = 2*W[14,k]+2*W[15,k]+2*W[16,k]
        U[10,k] = 1/2*U[9,k]
    return U, W

[U, W] = recurrence(U,W,x[:,0],k,mu_cb)

def propagate_int(x,k_max,mu_cb,t,t0,i):
    u0 = compute_u0(mu_cb, x[:,i-1])
    U = np.zeros((11,k_max+1))
    U[:,0] = u0 
    W = np.zeros((18,k_max+1))
    W[:,0] = compute_w0(U,x[:,i-1])
    xtest = np.zeros((11,1))
    [U, W] = recurrence(U,W,x[:,i-1],k_max,mu_cb)
    term = np.zeros((11,1))
    for k in range(1,k_max+1):
        term[:,0] = term[:,0] + 1/k*U[:,k-1]*(t-t0)**k
    xtest[:,0] = term[:,0] + x[:,i-1]
    return xtest

def propagate(U,W,x,k_max,mu_cb,t,steps):   
    for i in range(1,steps+1):
        t0_r = t[i-1]
        t_r = t[i] 
        xnew = propagate_int(x, k_max, mu_cb, t_r, t0_r,i)
        x[:,i] = xnew[:,0]
    return x

testing = propagate(U, W, x, k, mu_cb, t, steps)


#%%

plt.figure(0)
plt.plot(t,x[3], label ="$V_x$")
plt.plot(t,x[4], label="$V_y$")
plt.legend()
plt.show()


plt.figure(1)
plt.plot(t,x[0], label ="$R_x$")
plt.plot(t,x[1], label="$R_y$")
plt.plot(t,x[2], label="$R_z$")
plt.legend()
plt.show()

rmod = np.sqrt(x[7,:])

plt.figure(2)
plt.plot(t,rmod, label ="$radius$")
plt.legend()
plt.show()

#%%

print("Simulation Complete!")
