import numpy as np
import matplotlib.pyplot as plt

# Problem 12.4.7
# Time 
t0      = 0
tend    = 1
dt      = .2
t       = np.arange(t0, tend + dt, dt)

# Position
L       = 1
k       = 0.01
dx      = 0.01
x       = np.arange(0, L + dx, dx)

# Deflection
u       = np.zeros((len(x),len(t)))

for j in range(len(t)):
    for i in range(len(x)):
        u[i,j] = k*np.sin(2*np.pi*x[i])*np.cos(2*np.pi*t[j])
#     plt.figure()
#     plt.plot(x, u[:,j]) 
#     plt.ylabel('u(t,x)')
#     plt.xlabel('x')
#     plt.title(f't = {t[j]}')
#     plt.grid(True)
# plt.show()

plt.cycler(linestyle=['-', '-', '--', '-', '--', '-',],
                    color=['black', 'brown', 'brown', 'orange', 'orange', 'green']),
    
plt.figure()
plt.plot(x, u[:,0], label = "t = 0", linestyle="-") 
plt.plot(x, u[:,1], label = "t = 0.2", linestyle="--") 
plt.plot(x, u[:,2], label = "t = 0.4", linestyle="-.") 
plt.plot(x, u[:,3], label = "t = 0.6", linestyle=":") 
plt.plot(x, u[:,4], label = "t = 0.8", linestyle="-.") 
plt.plot(x, u[:,5], label = "t = 1.0", linestyle="--") 
plt.ylabel('u(t,x)')
plt.xlabel('x')
plt.title('Problem 12.4.7')
plt.grid(True)
plt.legend()
plt.show()

