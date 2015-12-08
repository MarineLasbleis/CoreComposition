#!/usr/bin/python

#### Composition of the core
### From Alfe 2002, Labrosse 2014, etc.



import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opti

np.random.seed() # initialize the random number generator with the current system time.

u = 1.660538921e-27
kb = 1.3806488e-23
eV = 1.60217656e-19
Tliq = 5500. #en Kelvin

# density of pure iron at ICB (Alfe 2002, EPSL)
rho_Fe_liq_ICB = 12968.
rho_Fe_sol_ICB = rho_Fe_liq_ICB*(1.017)


# partial atomic volume of elements (nu)  and  chemical potentials of impurities (from Alfe 2002 EPSL)
Fe = {'m' : 55.845*u,'nu_solid': 55.845*u/rho_Fe_sol_ICB, 'nu_liquid': 55.845*u/rho_Fe_liq_ICB}
S = {'m' : 32.06*u, 'nu_liquid': 6.65e-30, 'nu_solid': Fe["nu_solid"], 'Delta_mu' : 0.25*eV, 'lambda_l' : 6.15*eV, 'lambda_s' : 5.9*eV}
Si = {'m' : 28.085*u, 'nu_liquid': 6.65e-30, 'nu_solid': Fe["nu_solid"], 'Delta_mu' : 0.05*eV, 'lambda_l' : 3.6*eV, 'lambda_s' : 2.7*eV}
O = {'m' : 15.999*u,'nu_liquid': 4.25e-30, 'nu_solid': 4.65e-30, 'Delta_mu' : 2.6*eV, 'lambda_l' : 3.25*eV, 'lambda_s' : 0}



def Sum_mu(cs,cl,Delta_mu,lambdal,lambdas,T=Tliq):
    return Delta_mu+kb* T * np.log(cs/cl)-lambdal*cl+lambdas*cs


### Very simple model
## partition coeff = \Delta mu
## constant composition in IC
## Only Oxygene and Silicium

ccnO_OC = np.random.uniform(0.01,0.20,100)
ccnSi_OC = np.random.uniform(0.01,0.20,100)
ccnS_OC = np.random.uniform(0.01,0.20,100)

rho_OC_O = ((1-ccnO_OC)*Fe["m"]+ccnO_OC*O["m"]) / ((1-ccnO_OC)*Fe["nu_liquid"]+ccnO_OC*O["nu_liquid"])
rho_OC_Si = ((1-ccnSi_OC)*Fe["m"]+ccnSi_OC*Si["m"]) / ((1-ccnSi_OC)*Fe["nu_liquid"]+ccnSi_OC*Si["nu_liquid"])
rho_OC_S = ((1-ccnS_OC)*Fe["m"]+ccnS_OC*S["m"]) / ((1-ccnS_OC)*Fe["nu_liquid"]+ccnS_OC*S["nu_liquid"])



ccnO_IC = np.zeros(100)

for i, cl in enumerate(ccnO_OC):
 #   print cl, Sum_mu(1.e6,cl,O["Delta_mu"],O["lambda_l"],O["lambda_s"]), Sum_mu(0.3,cl,O["Delta_mu"],O["lambda_l"],O["lambda_s"])
    c = opti.brentq(Sum_mu, 1.e-10, 0.3 ,args=(cl,O["Delta_mu"],O["lambda_l"],O["lambda_s"]))
    ccnO_IC[i] = c
 #   print c, Sum_mu(c,cl,O["Delta_mu"],O["lambda_l"],O["lambda_s"])


ccnSi_IC = np.zeros(100)

for i, cl in enumerate(ccnSi_OC):
    c = opti.brentq(Sum_mu, 1.e-3, 0.3 ,args=(cl,Si["Delta_mu"],Si["lambda_l"],Si["lambda_s"]))
    ccnSi_IC[i] = c
    print c-cl, Sum_mu(c,cl,Si["Delta_mu"],Si["lambda_l"],Si["lambda_s"])

ccnS_IC=np.zeros(100)
for i, cl in enumerate(ccnS_OC):
    c=opti.brentq(Sum_mu, 1.e-3, 0.3 ,args=(cl,S["Delta_mu"],S["lambda_l"],S["lambda_s"]))
    ccnS_IC[i]=c


rho_IC_O = ((1-ccnO_IC)*Fe["m"]+ccnO_IC*O["m"]) / ((1-ccnO_IC)*Fe["nu_solid"]+ccnO_IC*O["nu_solid"])
rho_IC_Si = ((1-ccnSi_IC)*Fe["m"]+ccnSi_IC*Si["m"]) / ((1-ccnSi_IC)*Fe["nu_solid"]+ccnSi_IC*Si["nu_solid"])
rho_IC_S = ((1-ccnS_IC)*Fe["m"]+ccnS_IC*S["m"]) / ((1-ccnS_IC)*Fe["nu_solid"]+ccnS_IC*S["nu_solid"])


f, axarr = plt.subplots(3,1,sharex=True)
axarr[0].scatter(ccnO_OC, rho_OC_O, color='red')
axarr[0].scatter(ccnSi_OC, rho_OC_Si, color='blue')
axarr[0].scatter(ccnS_OC, rho_OC_S, color='green')
axarr[0].scatter(ccnO_OC, rho_IC_O, color='red')
axarr[0].scatter(ccnSi_OC, rho_IC_Si, color='blue')
axarr[0].scatter(ccnS_OC, rho_IC_S, color='green')
axarr[0].plot(np.linspace(0,0.20,3), 12166*np.ones(3))
axarr[0].set_ylabel('density (solid and liquid)')
axarr[1].scatter(ccnO_OC, ccnO_IC, color='red')
axarr[1].scatter(ccnSi_OC, ccnSi_IC, color='blue')
axarr[1].scatter(ccnS_OC, ccnS_IC, color='green')
axarr[1].plot(np.linspace(0,0.20,3), np.linspace(0,0.20,3))
axarr[1].set_ylabel('ccn in the solid')
#axarr[1].scatter(ccnO_OC, ccnO_OC* np.exp ((O["Delta_mu"]+O["lambda_l"]*ccnO_OC)/kb/Tliq))
axarr[2].scatter(ccnO_OC,(rho_IC_O-rho_OC_O)/rho_OC_O,color='red')
axarr[2].scatter(ccnSi_OC,(rho_IC_Si-rho_OC_Si)/rho_OC_Si,color='blue')
axarr[2].scatter(ccnS_OC,(rho_IC_S-rho_OC_S)/rho_OC_S,color='green')
axarr[2].plot(np.linspace(0,0.20,3),0.045*np.ones(3))
axarr[2].plot(np.linspace(0,0.20,3),0.06*np.ones(3))
axarr[2].set_ylabel('delta rho/rho_liq')
axarr[2].set_xlabel('ccn in liquid')
axarr[0].set_xlim([0,0.2])
axarr[0].set_ylim([11000,14000])
axarr[1].set_ylim([0, 0.2])
axarr[0].set_title('for O (red), Si (blue) and S (green)')


massO_OC=ccnO_OC*O["m"]/(ccnO_OC*O["m"]+(1-ccnO_OC)*Fe["m"])
massSi_OC=ccnSi_OC*Si["m"]/(ccnSi_OC*Si["m"]+(1-ccnSi_OC)*Fe["m"])
massS_OC=ccnS_OC*S["m"]/(ccnS_OC*S["m"]+(1-ccnS_OC)*Fe["m"])

f2,axarr2 =plt.subplots(3, 1, sharex=True)
axarr2[0].scatter(massO_OC, rho_OC_O, color='red')
axarr2[0].scatter(massSi_OC, rho_OC_Si, color='blue')
axarr2[0].scatter(massS_OC, rho_OC_S, color='green')
axarr2[0].scatter(massO_OC, rho_IC_O, color='red')
axarr2[0].scatter(massSi_OC, rho_IC_Si, color='blue')
axarr2[0].scatter(massS_OC, rho_IC_S, color='green')
axarr2[0].plot(np.linspace(0, 0.20, 3), 12166*np.ones(3))
axarr2[0].set_ylabel('density (solid and liquid)')
axarr2[1].scatter(massO_OC, ccnO_IC, color='red')
axarr2[1].scatter(massSi_OC, ccnSi_IC, color='blue')
axarr2[1].scatter(massS_OC, ccnS_IC, color='green')
axarr2[1].plot(np.linspace(0, 0.13, 3), np.linspace(0, 0.13, 3))
axarr2[1].set_ylabel('ccn in the solid')
#axarr[1].scatter(ccnO_OC, ccnO_OC* np.exp ((O["Delta_mu"]+O["lambda_l"]*ccnO_OC)/kb/Tliq))
axarr2[2].scatter(massO_OC, (rho_IC_O-rho_OC_O)/rho_OC_O, color='red')
axarr2[2].scatter(massSi_OC, (rho_IC_Si-rho_OC_Si)/rho_OC_Si, color='blue')
axarr2[2].scatter(massS_OC, (rho_IC_S-rho_OC_S)/rho_OC_S, color='green')
axarr2[2].plot(np.linspace(0, 0.20, 3), 0.045*np.ones(3))
axarr2[2].plot(np.linspace(0, 0.20, 3), 0.06*np.ones(3))
axarr2[2].set_ylabel('delta rho/rho_liq')
axarr2[2].set_xlabel('mass fraction in liquid')
axarr2[0].set_xlim([0, 0.13])
axarr2[0].set_ylim([11000, 14000])
axarr2[1].set_ylim([0, 0.13])
axarr2[0].set_title('for O (red), Si (blue) and S (green)')
plt.show()

