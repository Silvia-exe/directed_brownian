import numpy as np
import random as rnd
import matplotlib.pyplot as plt

#2D directed random walk of N particles

#Non-directed brownian motion, no boundary conditions, as control.
def true_brownian(v1):
    r0b = np.array([(0,0) for _ in range(N)])
    r_v = np.array([r0b])
    T_vb = np.arange(0, T,step = dt)
    pos_vb = np.array([r0b]) #Positions of particles

    r_oldb = r0b
    for i in range(0,T):
        r_oldb +=  np.array([r + v1*rnd.uniform(-1,1)*dt for r in r_oldb])
        r_v = np.append(r_v, r_oldb, axis = 0)
        true_brownian.append(r_old)
    for p in range(N):
        plt.plot(pos_vb[:, p, 0], pos_vb[:, p, 1])
    plt.xlabel("x-axis [x]")
    plt.ylabel("y-axis [y]")
    plt.title("2D true Brownian motion of "+ str(N) + " particles")
    plt.show()

#Implements reflectibe boundary conditions.
def reflective_bc(delta):
    if v0 == 0:
        return 0

    diff = abs(delta)//v0
    delta += -1*np.sign(delta)*abs(diff)*v0 #Periodic boundary conditions
    #delta *= 1-2*(abs(diff)%2) #Reflection?
    return delta

#Updates the position of the particle
def update_r(r_old, del_old):
    #print("position before update: ", r_old)
    deli = np.array([reflective_bc(d + v1*rnd.uniform(-1,1)) for d in del_old])
    ri = np.array([r_old[i] + deli[i]*dt for i in range(len(del_old))])
    #print("position after update: ", ri)
    return (ri, deli)

r0_0 = 0

T  = 5
dt = 0.1
N = 2
Lx = 10

v0 = 10.0 #Directed velocity, will describe self-propelled motion
tau_v = [0.1]

for tau in tau_v:
    v1 = (v0*dt)/tau #Random drift velocity, amplitude of random displacement

    r0 = np.random.uniform(-1,1,(N,2)) * Lx #Random initial conditions
    #r0 = np.array([(0,0) for _ in range(N)]) #Initial position as 0
    del0 = reflective_bc(np.random.uniform(-1,1,(N,2))*v0)
    #del0= np.array([(reflective_bc(rnd.uniform(-1,1)*v0), reflective_bc(rnd.uniform(-1,1)*v0)) for _ in range(N)])

    T_v = np.arange(0, T,step = dt)
    r_v = np.array([r0]) #Positions of particles
    del_v = np.array([del0]) #Ramdon fluctuations of particles
    #true_brownian = [r0] #True random brownian motion (DEBUG)
    posdiff_v = np.array([r0]) #For MSD WIP

    r_old = r0
    del_old = del0

    if N > 1:
        for t in range(len(T_v)-1):
            r_old, del_old = update_r(r_old, del_old)
            r_v = np.append(r_v, [r_old], axis = 0)
            del_v = np.append(del_v, [del_old], axis = 0)

            posdiff_v = np.append(posdiff_v, [(r0 - r_old)**2], axis = 0)

        msd_v = np.average(posdiff_v, axis = 1, keepdims= True) #MSD of particle ensemble
        plt.plot(msd_v, label = "$𝜏 = $ " + str(tau))

    else:
        for t in range(len(T_v)-1):
            r_old, del_old = update_r(r_old, del_old)
            r_v = np.append(r_v, [r_old], axis = 0)
            del_v = np.append(del_v, [del_old], axis = 0)


plt.xscale('log')
plt.yscale('log')
plt.xlabel("Time [t]")
plt.legend()
plt.ylabel("MSD [$x^2$] (Ensemble average)")
plt.title("Mean Square Displacement")
plt.show()

#pos_v[:, 0] #x and y of 0th particle
#pos_v[:, 1, 0] #x of pth particle'''

'''for p in range(N):
    plt.plot(r_v[:, p, 0], r_v[:, p, 1])
plt.xlabel("x-axis [x]")
plt.ylabel("y-axis [y]")
#plt.legend()
plt.title("2D directed Brownian motion of "+ str(N) + " particles")
plt.show()'''


'''plt.plot(T_v, del_v[:, 0], marker='.')
plt.xlabel("Time [t]")
plt.ylabel("$\delta(t)$ [x]")
plt.title("$\delta(t)$")
plt.show()'''
