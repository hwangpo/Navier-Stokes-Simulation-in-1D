#!/usr/bin/python
# -*- coding: utf-8 -*-

#----------------------------------------------------------------------------------------------------------
"""
Author : Mohamed

Title : Implementation des equations de saint-venant

Description : Simulation des vagues dans un domaine [0,1] avec les equations de saint-venant

"""
#-----------------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

#Params
g = 9.8

def F(U):
	"""
	U[0] : hauteur
	U[1] : debit
	"""

	return np.array([U[1],(U[1]**2)/U[0] + g*(U[0]**2)/2])

def f(U_g,U_d):
	"""
	U_g/d : vecteur a gauche/droite
	"""

	#vitesse
	v_g = U_g[1]/U_g[0]
	v_d = U_d[1]/U_d[0]

	c = np.amax([np.abs(v_g)+np.sqrt(g*U_g[0]), np.abs(v_d)+np.sqrt(g*U_d[0])])

	return (F(U_g)+F(U_d))/2 - (c/2)*(U_d - U_g)



#----------------------------------------------------------------------------------------------------------
# Mayage
#-------------------------

N = 50
dx = 1.0/N

x = np.linspace(0, 1, num=N+1)  # le nombre d'element est N+1 car on a N intervalle

#[End]----------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------
# Conditions initiales et de bords
#-------------------------

v_0 = np.zeros(N+1) #vitesse nulle a t = 0


#Condition discontinue sur les hauteurs
# h_0 = np.ones(N+1)
# h_0[N//2:] = np.ones(len(h_0[N//2:])) * 20

#Condition continue
h_0 = np.linspace(1,20,N+1)

U_0 = np.array([ h_0, h_0*v_0 ])

#[End]----------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
# Definition des tableaux
#-------------------------

U_tab = np.array([U_0]) #tableau de vecteur qui va contenir les vecteurs U (qui contient les donnees de tout le maillage) a tout instant.

time = np.zeros(1) #Tableau des instants

#[End]----------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
# Implementation
#-------------------------

def lambda1(v,h):
	return v - np.sqrt(g*h)

def lambda2(v,h):
	return v + np.sqrt(g*h)	

#debut de l'evolution
n = 0 #--> n designe le temps

while n < 1000:

	h_n = U_tab[n][0] #tableau hauteurs h_n
	q_n = U_tab[n][1] #tableau debits q_n

	v_n = q_n/h_n #tableau vitesses v_n

	m = 2 * np.max([np.abs(lambda1(v_n,h_n)), np.abs(lambda2(v_n,h_n))]) 
	dt = dx/m

	time = np.append(time, time[n]+dt) #Ajout dans le tableau de temps

	# U_nplus1 = U_tab[n]  #--> fait gaffe la, car U_nplus1 prend la meme adress que U_tab[n] !!!!!
	U_nplus1 = np.array([np.ones(N+1), np.ones(N+1)]) #Tableau qui contient les valeurs de vecteurs U a l'instant n+1 sur tout le maillage

	for i in range(1,N): # i -> designe l'endroit sur le maillage

	#Approche: Vitesse nulle au bord (pour faire une reflexion)

		u_in = np.array([h_n[i], q_n[i]]) #u_in vecteur a l'endroit i qui contient les donnees h_n[i],q_n[i] a l'instant n et a l'endroit i

		u_ig = np.array([h_n[(i-1)], q_n[(i-1)]]) #definition de u_gauche a l'instant n avec condition periodique
		u_id = np.array([h_n[(i+1)], q_n[(i+1)]])

		u_i_nplus1 = u_in - (dt/dx)*(f(u_in, u_id)-f(u_ig, u_in)) #def de vecteur U a l'instant n+1 a l'endroit i

		#Stocker les params dans le nouveau tableau de vecteur a l'instant n+1 a l'endroit i
		U_nplus1[0][i] = u_i_nplus1[0] 
		U_nplus1[1][i] = u_i_nplus1[1]

	#i=0

	U_nplus1[0][0] = h_n[1]
	U_nplus1[1][0] = -0.7*q_n[1] #Vitesse avec un coefficient de restitution

	#i=N

	U_nplus1[0][N] = h_n[N-1] 
	U_nplus1[1][N] = -0.7*q_n[N-1] #Vitesse avec un coefficient de restitution




	U_tab = np.concatenate((U_tab, [U_nplus1])) #Ajout des nouveaux donnees de l'instant n+1

	n = n+1 #incrementation de temps


#[End]----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
# Animation de plot
#-------------------------


# plt.ion() # mode interactif
# plt.figure(1) # nouvelle figure    

# plt.plot(x, U_tab[0][0])
# plt.show()

# plt.pause(0.01) # nécessaire (bug matplotlib)
# plt.hold(False) # Ne pas superposer les figures pour animer le mouvement

# for j in range(1,n):
   
#     plt.plot(x, U_tab[j][0])  # on redessine
#     plt.ylim([0, 21]) #limits sur l'axis y
#     plt.pause(0.00001) # on attend un peu

# plt.ioff()        # fin du mode interactif
# plt.show()        #nécessaire

import matplotlib.animation as animation  #Tutorial http://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

fig = plt.figure()

plt.ylim([-80, 80]) #limits sur l'axis y

line1, = plt.plot(x, U_tab[0][0], lw=2)
line2, = plt.plot(x, U_tab[0][1], lw=2)

def animate(j):
    y1 = U_tab[j][0]
    x1 = x
    line1.set_data(x1, y1)

    y1 = U_tab[j][1]
    line2.set_data(x1, y1)

    return line1,line2
        
def init():
    line1.set_data(x, U_tab[0][0])
    line2.set_data(x, U_tab[0][1])
    return line1,line2
          
ani = animation.FuncAnimation(fig, animate, init_func = init, frames = n, interval = 5, blit = False)  

ani.save('wave_motion.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()

#[End]----------------------------------------------------------------------------------------------------