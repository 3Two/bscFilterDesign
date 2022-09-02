#!/usr/bin/env python

import scipy.interpolate as interpolate
import interp_funcs as interp_funcs
import numpy as np

# Attenuation = 1- absorption
def attenuation(E,element,density,thickness): # E in MeV, element as string, density in g/cc, thickness in mm

    attenuation_data_path = './Attenuation_Data/'
    thickness            *= 0.1 # convert mm to cm
    mass_thickness        = thickness * density # mass tickness in g/cm2
    attenData             = np.loadtxt(attenuation_data_path+element+'.dat',comments="#") # Load in attenuation data
    photon_MeV, mu_rho    = attenData[:,0], attenData[:,1]

    filename = element+'_interp.dat'
    file = open(attenuation_data_path+filename, 'r')
    interp_edge = []
    interp_kind = []
    for line in file:
        p = line.split()
        interp_edge.append(int(p[0]))
        interp_kind.append(p[1])
    file.close()

    if(interp_edge[0]==0):
        interp = interp_funcs.no_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind)
    elif(interp_edge[0]==1):
        interp = interp_funcs.one_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind)
    elif(interp_edge[0]==2):
        interp = interp_funcs.two_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind)
    elif(interp_edge[0]==4):
        interp = interp_funcs.four_edge(E,photon_MeV,mu_rho,interp_edge,interp_kind)

    attenuation_factor = np.exp(-(interp*mass_thickness))

    """Plotting Check"""
    if(0):
        import matplotlib.pyplot as plt
        plt.title(element)
        plt.plot(photon_MeV*1e3, mu_rho,label='NIST')
        plt.plot(E*1e3,interp,label='Interp')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Photon Energy (keV)')
        plt.ylabel(r'$\mu /\rho$')
        plt.legend(loc='best')
        plt.show()

    #plotting check
    if(0):
        import matplotlib.pyplot as plt
        plt.plot(E*1e3,attenuation_factor,label=element)
        plt.legend(loc='best')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Photon Energy (keV)')
        plt.ylabel(r'$\mu /\rho$')
        plt.show()

    return attenuation_factor
