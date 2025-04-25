import numpy as np
import random as rnd
import matplotlib.pyplot as plt

#1D directed random walk of N particles

#Implements reflectibe boundary conditions.
def reflective_bc(delta):
    if v0 == 0:
        return delta

    diff = abs(delta)//v0
    delta += -1*np.sign(delta)*abs(diff)*v0 #Periodic boundary conditions
    #delta *= 1-2*(abs(diff)%2) #Reflection?
    return delta

#Updates the position of the particle
def update_r(r_old, del_old):
    deli = np.array([reflective_bc(d + v1*rnd.uniform(-1,1)) for d in del_old])
    ri = np.array([r_old[i] + deli[i]*dt for i in range(len(del_old))])
    return (ri, deli)


T = 10**3
dt = 0.5
N = 5

v0 = 0.0 #Directed velocity, will describe self-propelled motion
tau_v = [100]
i = 0

for tau in tau_v:
    print("Running for tau = ", tau)
    print(str(100*i/len(tau_v)), "% done." )
    v1 = (v0*dt)/tau #Random drift velocity, amplitude of random displacement

    r0 = np.zeros(N)
    del0= np.array([reflective_bc(rnd.uniform(-1,1)*v0) for _ in range(N)])

    T_v = np.arange(0, T,step = dt)
    pos_v = np.array([r0]) #Positions of particles
    del_v = np.array([del0]) #Ramdon fluctuations of particles
    posdiff_v = np.array([r0]) #For MSD WIP

    r_old = r0
    del_old = del0

    if N > 1:
        for t in range(len(T_v)-1):
            r_old, del_old = update_r(r_old, del_old)
            pos_v = np.append(pos_v, [r_old], axis = 0)
            del_v = np.append(del_v, [del_old], axis = 0)
            posdiff_v = np.append(posdiff_v, [(r0 - r_old)**2], axis = 0)

        msd_v = np.average(posdiff_v, axis = 1, keepdims= True) #MSD of particle ensemble
        plt.plot(msd_v, label = "$𝜏 = $ " + str(tau))

    else:
        for t in range(len(T_v)-1):
            r_old, del_old = update_r(r_old, del_old)
            pos_v = np.append(pos_v, [r_old], axis = 0)
            del_v = np.append(del_v, [del_old], axis = 0)
    i += 1

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Time [t]")
plt.legend()
plt.ylabel("MSD [$x^2$] (Ensemble average)")
plt.title("Mean Square Displacement")
plt.show()

'''r_old = r0
for i in range(0,T):
    r_old += rnd.uniform(-1,1)*v1
    true_brownian.append(r_old)
plt.plot(true_brownian, marker='.')
plt.title("true_brownian")
plt.xlabel("Time [t]")
plt.ylabel("Position 1D [x]")
plt.show()'''

sel_p = np.array([rnd.randint(0,N) for _ in range(100)])
for p in sel_p:
    plt.plot(T_v, pos_v[:, p])
#plt.axline((0,r0), slope = del0, color = "red", lw = 1, label ="Direction of dir. movement.")
plt.xlabel("Time [t]")
plt.ylabel("Position 1D [x]")
#plt.legend()
plt.title("1D directed Brownian motion of "+ str(N) + " particles")
plt.show()

'''plt.plot(T_v, del_v[:, 0], marker='.')
plt.xlabel("Time [t]")
plt.ylabel("$\delta(t)$ [x]")
plt.title("$\delta(t)$")
plt.show()'''
