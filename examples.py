# -*- coding: utf-8 -*-
"""
Examples of plots and calculations using the tmm package.
"""

#make division of integers work as expected
from __future__ import division

from tmm import *
from numpy import pi, linspace
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# "5 * degree" is 5 degrees expressed in radians
# "1.2 / degree" is 1.2 radians expressed in degrees
degree = pi/180

def sample1():
    """
    Here's a thin non-absorbing layer, on top of a thick absorbing layer, with
    air on both sides. Plotting reflected intensity versus wavenumber, at two
    different incident angles.
    """
    d_list = [inf,100,300,inf] #in nm
    n_list = [1,2.2,3.3+0.3j,1]
    plotpoints=400
    ks=linspace(0.0001,.01,plotpoints) #in nm^-1
    Rnorm=zeros(plotpoints)
    R45=zeros(plotpoints)
    for i in range(plotpoints):
        norm_data = coh_tmm('s',n_list, d_list, 0, 1/ks[i]) #in nm
        s45_data = coh_tmm('s',n_list, d_list, 45*degree, 1/ks[i]) #in nm
        p45_data = coh_tmm('p',n_list, d_list, 45*degree, 1/ks[i]) #in nm
        Rnorm[i]=norm_data['R']
        R45[i]=(s45_data['R']+p45_data['R'])/2
    kcm = ks * 1e7 #ks in cm^-1 rather than nm^-1
    plt.figure()
    plt.plot(kcm,Rnorm,'blue',kcm,R45,'purple')
    plt.xlabel('k (cm$^{-1}$)')
    plt.ylabel('Fraction reflected')
    plt.title('Reflection of unpolarized light at 0$^\circ$ incidence (blue), '
                '45$^\circ$ (purple)')

def sample2():
    """
    Here's the transmitted intensity versus wavelength through a single-layer
    film which has some complicated wavelength-dependent index of refraction.
    (I made these numbers up, but in real life they could be read out of a
    graph / table published in the literature.) Air is on both sides of the
    film, and the light is normally incident.
    """
    #index of refraction of my material: wavelength in nm versus index.
    material_nk_data = array([[200, 2.1+0.1j],
                              [300, 2.4+0.3j],
                              [400, 2.3+0.4j],
                              [500, 2.2+0.4j],
                              [750, 2.2+0.5j]])
    material_nk_fn = interp1d(material_nk_data[:,0].real,
                              material_nk_data[:,1], kind='quadratic')
    d_list = [inf,300,inf] #in nm
    lambda_list = linspace(200,750,400) #in nm
    T_list = []
    for lambda_vac in lambda_list:
        n_list = [1, material_nk_fn(lambda_vac), 1]
        T_new = coh_tmm('s',n_list,d_list,0,lambda_vac)['T']
        T_list.append(T_new)
    plt.figure()
    plt.plot(lambda_list,T_list)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Fraction of power transmitted')
    plt.title('Transmission at normal incidence')
        

def sample3():
    """
    Here is a calculation of the psi and Delta parameters measured in
    ellipsometry. This reproduces Fig. 1.14 in Handbook of Ellipsometry by
    Tompkins, 2005.
    """
    n_list=[1,1.46,3.87+0.02j]
    plotpoints=100
    ds=linspace(0,1000,plotpoints) #in nm
    psis=zeros(plotpoints)
    Deltas=zeros(plotpoints)
    for i in range(plotpoints):
        e_data=ellips(n_list, [inf,ds[i],inf], 70*degree, 633) #in nm
        psis[i]=e_data['psi']/degree # angle in degrees
        Deltas[i]=e_data['Delta']/degree # angle in degrees
    plt.figure()
    plt.plot(ds,psis,ds,Deltas)
    plt.xlabel('SiO2 thickness (nm)')
    plt.ylabel('Ellipsometric angles (degrees)')
    plt.title('Ellipsometric parameters for air/SiO2/Si, varying '
            'SiO2 thickness.\n' 
            '@ 70$^\circ$, 633nm. '
            'Should agree with Handbook of Ellipsometry Fig. 1.14')

def sample4():
    d_list = [inf, 100, 300, inf] #in nm
    n_list = [1, 2.2+0.2j, 3.3+0.3j, 1]
    th_0=pi/4
    lam_vac=400
    pol='p'
    coh_tmm_data = coh_tmm(pol,n_list,d_list,th_0,lam_vac)
    
    plotpoints=1000
    ds = linspace(0,400,plotpoints) #position in structure
    poyn=zeros(plotpoints)
    absor=zeros(plotpoints)
    for i in range(plotpoints):
        layer, d_in_layer = find_in_structure_with_inf(d_list,ds[i])
        data=position_resolved(layer,d_in_layer,coh_tmm_data)
        poyn[i] = data['poyn']
        absor[i] = data['absor']
    plt.figure()
    plt.plot(ds,poyn,'blue',ds,200*absor,'purple')
    plt.xlabel('depth (nm)')
    plt.ylabel('AU')
    plt.title('Local absorption (purple), Poynting vector (blue)')

