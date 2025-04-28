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

#Acceptance function for movement
#Will stop a particle if it would intersect another
def p_acceptance(r_old, r_new):
    for p in range(N):
        for k in range(N):
            if (np.sqrt((r_new[p,0] - r_old[k,0])**2 + (r_new[p,1] - r_old[k,1])**2)) < 2*sigma:
                if k != p:
                    print("old movement ", r_old[p])
                    print("new movement ", r_new[p])
                    print("could collide with  ", r_old[k])
                    r_new[p] = r_old[p]
    return r_new

#Updates the position of the particle
def update_r(r_old, del_old):
    #print("position before update: ", r_old)
    del_new = np.array([reflective_bc(d + v1*rnd.uniform(-1,1)) for d in del_old])
    r_new = np.array([r_old[i] + del_new[i]*dt for i in range(len(del_old))])
    #print("position after update: ", ri)
    return (r_new, del_new)

def initialize_r0():
    


#Simulation variables
T  = 1000 #Total simulation time
dt = 0.1 #dt for T
N = 50 #Number of particles
sigma = 1 #Particle radius and system length scale
Lx = 50*sigma #Box size

v0 = 1.0 #Directed velocity, will describe self-propelled motion
tau_v = [2.1] #Persistance time

#Figure settings
fig, ax = plt.subplots()
ax.set_xlim(0, Lx)
ax.set_ylim(0, Lx)
plt.gca().set_aspect('equal', adjustable='box')

#Position update loop for each specified tau
for tau in tau_v:
    v1 = (v0*dt)/np.sqrt(tau) #Random drift velocity, amplitude of random displacement

    r0 = np.random.rand(N,2) * Lx #Random initial conditions CHECK THEY ARENT OVERLAPPING!!

    del0 = reflective_bc(np.random.uniform(-1,1,(N,2))*v0) #Random initial delta

    T_v = np.arange(0, T,step = dt) #Time vector

    r_v = np.array([r0]) #Positions of particles
    del_v = np.array([del0]) #Random fluctuations of particles

    #posdiff_v = np.array([r0]) #For MSD WIP

    r_old = r0
    del_old = del0

    for p in r0:
        circles = ax.add_artist(plt.Circle((p[0], p[1]), sigma, fill= False))

    plt.show()
    #scatter = ax.scatter(r0[:, 0], r0[:, 1], marker='o', s=np.pi*sigma**2)

    if N > 1:
        for t in range(len(T_v)-1):
            r_new, del_new = update_r(r_old, del_old)

            r_new = p_acceptance(r_old, r_new)

            #Detect boundary crossing and implement PBC
            for p in range(N):
                for i in range(2):
                    if r_new[p,i] < 0:
                        r_new[p,i] += Lx
                    elif r_new[p,i] > Lx:
                        r_new[p,i] -= Lx

            #Saving of particle positions
            r_v = np.append(r_v, [r_new], axis = 0)
            del_v = np.append(del_v, [del_new], axis = 0)

            #posdiff_v = np.append(posdiff_v, [(r0 - r_old)**2], axis = 0)

            #Plotting positions as animation (fun)
            circles.set_center(r_new)
            plt.pause(0.001)

            r_old = r_new
            del_old = del_new

        plt.show()
        #msd_v = np.average(posdiff_v, axis = 1, keepdims= True) #MSD of particle ensemble
        #plt.plot(msd_v, label = "$𝜏 = $ " + str(tau))

    else:
        for t in range(len(T_v)-1):
            r_old, del_old = update_r(r_old, del_old)
            r_v = np.append(r_v, [r_old], axis = 0)
            del_v = np.append(del_v, [del_old], axis = 0)


'''plt.xscale('log')
plt.yscale('log')
plt.xlabel("Time [t]")
plt.legend()
plt.ylabel("MSD [$x^2$] (Ensemble average)")
plt.title("Mean Square Displacement")
plt.show()'''

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
