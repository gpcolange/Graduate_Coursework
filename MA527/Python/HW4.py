import numpy as np
import matplotlib.pyplot as plt

def heaviside(t,a):
    # Heaviside Step Function u(t - a)
    if t < a:
        u = 0
    else:
        u = 1
        
    return u

# Problem 6.3.13
t0      = 0
tend    = 5
dt      = .005
t       = np.arange(t0, tend + dt, dt)
f       = np.zeros(len(t))

for i in range(len(t)):
    f[i] = 2*np.sin(3*t[i])*(1 +  heaviside(t[i],np.pi))

fig1 = plt.figure(1)
plt.plot(t, f) 
plt.ylabel('f(t)')
plt.xlabel('Time [s]')
plt.title('Problem 6.3.13')
plt.grid(True)

# Problem 6.3.16
f       = np.zeros(len(t))

for i in range(len(t)):
    f[i] = (1/2)*heaviside(t[i],1)*(np.exp(2*t[i] - 2) - np.exp(-2*t[i] + 2)) + \
            - (1/2)*heaviside(t[i],3)*(np.exp(2*t[i] - 6) - np.exp(-2*t[i] + 6))

fig2 = plt.figure(2)
plt.plot(t, f) 
plt.ylabel('f(t)')
plt.xlabel('Time [s]')
plt.title('Problem 6.3.16')
plt.grid(True)

# Problem 6.4.3
f       = np.zeros(len(t))

for i in range(len(t)):
    if t[i] > np.pi:
        f[i] = 8*np.cos(2*t[i]) + 0.5*np.sin(2*t[i])
    else: 
        f[i] = 8*np.cos(2*t[i])

fig3 = plt.figure(3)
plt.plot(t, f)
plt.ylabel('f(t)')
plt.xlabel('Time [s]')
plt.title('Problem 6.4.3')
plt.grid(True)

# Problem 6.4.10
f       = np.zeros(len(t))

for i in range(len(t)):
    if t[i] < np.pi/2:
        f[i] = 0
    elif t[i] > np.pi/2 and t[i] < np.pi:
        f[i] = np.exp(-2*(t[i] - np.pi/2)) - np.exp(-3*(t[i] - np.pi/2))
    else:
        f[i] = np.exp(-2*(t[i] - np.pi/2)) - np.exp(-3*(t[i] - np.pi/2)) - np.cos(t[i] - np.pi)/10 \
              - np.sin(t[i] - np.pi)/10 + (4/10)*np.exp(-2*(t[i] - np.pi)) - (3/10)*np.exp(-3*(t[i] - np.pi))

fig4 = plt.figure(4)
plt.plot(t, f)
plt.ylabel('f(t)')
plt.xlabel('Time [s]')
plt.title('Problem 6.4.10')
plt.grid(True)
plt.show()