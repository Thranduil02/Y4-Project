#Project Graphs for D-Bion

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets 
from scipy.optimize import fsolve
import sympy as sp

#============================================= CONSTANTS ================================================

beta=1
rho=1
M_p=1
r_0=1
Lambda=1

#slider limits
betamin = 0
betamax = 10
rhomin = 0
rhomax = 10
r_0min = 0
r_0max = 10
lambdamin = 0
lambdamax = 10

inital_slider_val = 1

r_v_spherical=np.sqrt((4*beta*(r_0**3)*rho)/(12*M_p*(Lambda**2)))
r_v_cylindrical=((r_0**2)*rho*beta)/(2*M_p*(Lambda**2))
r_star=(M_p*(Lambda**2))/(beta*rho)

#planar source

r = np.linspace(0,10,1000)

#=========================================== FUNCTIONS ==================================================

def phidash_planar(r, r_0, beta, rho, M_p):
    """
    Function to calculate phidash for planar source
    """
    phidash = np.zeros_like(r)
    # Initialize phidash as an array of zeros
    
    for idx, pos in enumerate(r):
        if pos < r_0:
            phidash[idx]=(Lambda**2)/(np.sqrt(1+((r_star**2)/(pos**2))))
        else:
            phidash[idx]=(Lambda**2)/(np.sqrt(1+((r_star**2)/(r_0**2))))
            
    return phidash

#spherical source

def phidash_spherical(r, r_0, beta, rho, M_p, Lambda, r_v_spherical):
    """
    Function to calculate phidash for spherical source
    """
    phidash = np.zeros_like(r)
    # Initialize phidash as an array of zeros
    
    for idx, pos in enumerate(r):
        if pos < r_0:
            phidash[idx] = (Lambda**2)/(np.sqrt(1+((r_0**6)/((pos**2)*(r_v_spherical**4)))))
        else:
            phidash[idx] = (Lambda**2)/(np.sqrt(1+((pos**4)/((r_v_spherical**4)))))
            
    return phidash

#cylindrical source

def phidash_cylindrical(r, r_0, beta, rho, M_p, Lambda, r_v_cylindrical):
    """
    Function to calculate phidash for spherical source
    """
    phidash = np.zeros_like(r)
    # Initialize phidash as an array of zeros
    
    for idx, pos in enumerate(r):
        if pos < r_0:
            phidash[idx] = (Lambda**2)/(np.sqrt(1+((r_0**4)/((pos**2)*(r_v_cylindrical**2)))))
        else:
            phidash[idx] = (Lambda**2)/(np.sqrt(1+((pos**2)/((r_v_cylindrical**2)))))
            
    return phidash


#============================================ PLOTTING ======================================================

def update(val):
    """
    Function to update the plot when the sliders are moved
    """
    beta = slider_beta.val
    rho = slider_rho.val
    r_0 = slider_r_0.val
    Lambda = slider_Lambda.val
    
    r_v_spherical=np.sqrt((4*beta*(r_0**3)*rho)/(12*M_p*(Lambda**2)))
    r_v_cylindrical=((r_0**2)*rho*beta)/(2*M_p*(Lambda**2))
    r_star=(M_p*(Lambda**2))/(beta*rho)
    
    axes["A"].clear()
    axes["A"].plot(r,phidash_planar(r, r_0, beta, rho, M_p), color='red')
    axes["A"].set_xlabel('$r$', fontsize=14)
    axes["A"].set_ylabel("$\phi'$", fontsize=14)
    axes["A"].set_title('Planar Source')
    axes["A"].legend()

    axes["B"].clear()
    axes["B"].plot(r, phidash_spherical(r, r_0, beta, rho, M_p, Lambda, r_v_spherical), color='blue')
    axes["B"].set_xlabel('$r$', fontsize=14)
    axes["B"].set_ylabel("$\phi'$", fontsize=14)
    axes["B"].set_title('Spherical Source')
    axes["B"].legend()

    axes["C"].clear()
    axes["C"].plot(r, phidash_cylindrical(r, r_0, beta, rho, M_p, Lambda, r_v_cylindrical), color='green')
    axes["C"].set_xlabel('$r$', fontsize=14)
    axes["C"].set_ylabel("$\phi'$", fontsize=14)
    axes["C"].set_title('Cylindrical Source')
    axes["C"].legend()

    axes["X"].clear()
    axes["X"].axis('off')
    axes["Z"].clear()
    axes["Z"].axis('off')
    axes["Z"].set_title('The D-Bion for Planar, Spherical and Cylindrical Sources')

    fig.canvas.draw()

# Initialise Plot
fig = plt.figure()

axes = fig.subplot_mosaic(
    """
    ZZZ
    ABC
    ABC
    XXX
    """
)

axes["A"].plot(r,phidash_planar(r, r_0, beta, rho, M_p), color='red')
axes["A"].set_xlabel('$r$', fontsize=14)
axes["A"].set_ylabel("$|\phi'|$", fontsize=14)
axes["A"].set_title('Planar Source')
axes["A"].legend()

axes["B"].plot(r, phidash_spherical(r, r_0, beta, rho, M_p, Lambda, r_v_spherical), color='blue')
axes["B"].set_xlabel('$r$', fontsize=14)
axes["B"].set_ylabel("$|\phi'|$", fontsize=14)
axes["B"].set_title('Spherical Source')
axes["B"].legend()

axes["C"].plot(r, phidash_cylindrical(r, r_0, beta, rho, M_p, Lambda, r_v_cylindrical), color='green')
axes["C"].set_xlabel('$r$', fontsize=14)
axes["C"].set_ylabel("$|\phi'|$", fontsize=14)
axes["C"].set_title('Cylindrical Source')
axes["C"].legend()

axes["X"].axis('off')
axes["Z"].axis('off')
axes["Z"].set_title('The D-Bion for Planar, Spherical and Cylindrical Sources')

#============================================ SLIDERS ================================================

#slider axes
slider_width = 0.01
slider_length = 0.65

beta_slider_ax = plt.axes([0.1, 0.20, slider_length, slider_width])
rho_slider_ax = plt.axes([0.1, 0.15, slider_length, slider_width])
r_0_slider_ax = plt.axes([0.1, 0.10, slider_length, slider_width])
Lambda_slider_ax = plt.axes([0.1, 0.05, slider_length, slider_width])

#sliders
slider_beta = widgets.Slider(beta_slider_ax, r'$\beta$', betamin, betamax, valinit=inital_slider_val)
slider_rho = widgets.Slider(rho_slider_ax, r'$\rho_0$', rhomin, rhomax, valinit=inital_slider_val)
slider_r_0 = widgets.Slider(r_0_slider_ax, r'$r_0$', r_0min, r_0max, valinit=inital_slider_val)
slider_Lambda = widgets.Slider(Lambda_slider_ax, r'$\Lambda$', lambdamin, lambdamax, valinit=inital_slider_val)

#slider callbcks
slider_beta.on_changed(update)
slider_rho.on_changed(update)
slider_r_0.on_changed(update)
slider_Lambda.on_changed(update)
                                
plt.show()

