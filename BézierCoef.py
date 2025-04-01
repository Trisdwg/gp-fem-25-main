import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

Tmin = np.pi/2 - 0.5
Tmax = 1.5707963
n = 10
p = 3
T = np.linspace(Tmin,Tmax+(p-1)*(Tmax-Tmin)/n,n+p)
tau = np.linspace(T[p],T[n-p],100)

Rin = 25

def x(t):
    return Rin*np.cos(t) - np.sin(t)*Rin*(np.pi/2 - t)

def y(t):
    return Rin*np.sin(t) + np.cos(t)*Rin*(np.pi/2 - t)
    

def Spline(t, i, p):
    if(p == 0):
        if(T[i] <= t and t < T[i+1]):
            return 1
        else:
            return 0
    else:
        return (t-T[i])/(T[i+p]-T[i])*Spline(t,i,p-1) + (T[i+1+p]-t)/(T[i+1+p]-T[i+1])*Spline(t,i+1,p-1)

A = np.empty((n-p,n-p))
b = np.empty(n-p)
c = np.empty(n-p)
for k in range(n-p):
    for l in range(n-p):
        A[k][l] = integrate.quad(lambda x: Spline(x,k,p)*Spline(x,l,p), T[p], T[n-p])[0]
    b[k] = integrate.quad(lambda t: Spline(t,k,p)*x(t), T[p], T[n-p])[0]
    c[k] = integrate.quad(lambda t: Spline(t,k,p)*y(t), T[p], T[n-p])[0]
    
X = np.linalg.solve(A,b)
Y = np.linalg.solve(A,c)
xh = np.zeros_like(tau)
resx = np.empty_like(tau)
for m in range(len(tau)):
    resx[m] = x(tau[m])
    for i in range(n-p):
        xh[m] += Spline(tau[m],i,p)*X[i]
yh = np.zeros_like(tau)
resy = np.empty_like(tau)
for m in range(len(tau)):
    resy[m] = y(tau[m])
    for i in range(n-p):
        yh[m] += Spline(tau[m],i,p)*Y[i]

fig, ax = plt.subplots()
#ax.set_aspect('equal', adjustable='box')
plt.plot(xh,yh)
plt.plot(resx,resy,'o')
plt.plot(X,Y,'o')
plt.grid(True)
plt.show()