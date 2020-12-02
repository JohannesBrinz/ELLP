''' Versuch ELTT
 Johannes Brinz & Caterina Vanelli
 Betreuer: Eduard Bykov
 Datum: 13.11.2020
 Ort: Rossendorf Forschungszentrum'''

import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
from scipy import optimize
import math
import matplotlib.font_manager as fm


# Cu_Si_1 Metall und Halbleiter
#Datenimport
Cu_Si_1 = pd.read_csv("Daten/Cu_and_Si_1.dat", sep = "\t", header = 4, \
    names = ["Zeit [s]","T [K]", "R_Probe_1 [Ohm]", "R_Thermom [Ohm]", "R_Probe_2 [Ohm]"])


#Plot
plt.errorbar(Cu_Si_1["T [K]"], Cu_Si_1["R_Probe_1 [Ohm]"], \
    xerr = 0, yerr = 0, c='black', ecolor='red', fmt='.')

plt.title('Cu_Si_1')
plt.xlabel('T in K')
plt.ylabel('R_Probe_1 in Ohm')

plt.savefig('Plots/Cu_and_Si_1.png', dpi=300)

plt.clf()



#Nb_and_Si: Superconductor
#Without Magnetic field
#Data import
Nb_Si_2 = pd.read_csv("Daten/Nb_and_Si_2.dat", sep = "\t", header = 4, \
    names = ["Zeit[s]","T[K]", "R_Probe_1[Ohm]", "R_Th[Ohm]", "R_Probe_2[Ohm]"])


#Fit Logistic Function
def Logis(T, G, k, theta, F):
    return G * 1 /( 1 + np.exp(-k * (T - theta)) ) + F

params, params_cov = optimize.curve_fit(Logis, Nb_Si_2["T[K]"][2760:2803], Nb_Si_2["R_Probe_1[Ohm]"][2760:2803], \
p0 = [0.03, 1, 9.5, 0.017])
perr = np.sqrt(np.diag(params_cov)) #Error of the fit parameters

print("G = ", params[0], "k = ", params[1], "F = ", params[3] )

#Plot
plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
plt.gcf().subplots_adjust(top=0.9)

plt.errorbar(Nb_Si_2["T[K]"], Nb_Si_2["R_Probe_1[Ohm]"], c='blue', fmt='x')
plt.errorbar(np.linspace(7,11, 1000), Logis(np.linspace(7,11, 1000), \
    params[0], params[1], params[2], params[3]), label='Fitted function', c= 'black')

plt.errorbar(params[2], Logis(params[2], params[0], params[1], params[2], params[3]), xerr = perr[2], capsize=3, c = "red", fmt = "+" )

plt.xlim(8, 10.5)
plt.ylim(0.015, 0.0325)

plt.title('Super Conductor without Magnetic Field', fontsize = 20)
plt.xlabel('T [K]', fontsize = 14)
plt.ylabel('$R_{SuperCond}$ [$\Omega$]', fontsize = 14)

plt.text( 8.25, 0.027, "$T_c$ = " + str(round(params[2], 2)) + \
    " $\pm$ " + str(round(perr[2], 3)), fontsize = 15 )
plt.legend(['Measured Resistivity',  "Fitted Plot", "$T_c$"], fontsize = 13)

plt.savefig('Plots/Nb_and_Si_2.png', dpi=300)
plt.clf()




#With Magnetic field
#Data import
Nb_Si_MF = pd.read_csv("Daten/Nb_and_Si_MF.dat", sep = "\t", header = 4, \
    names = ["Zeit[s]","T[K]", "R_Probe_1[Ohm]", "R_Th[Ohm]", "R_Probe_2[Ohm]"])

#Fit Logistic Function
#1 Ampere
params_1, params_cov_1 = optimize.curve_fit(Logis, Nb_Si_MF["T[K]"][780:920],\
    Nb_Si_MF["R_Probe_1[Ohm]"][780:920], p0 = [0.03, 1, 9, 0.017])
perr_1 = np.sqrt(np.diag(params_cov_1))

print("1 Ampere", params_1[2], "+-", perr_1[2])

#2 Ampere
params_2, params_cov_2 = optimize.curve_fit(Logis, Nb_Si_MF["T[K]"][1106:1271], \
    Nb_Si_MF["R_Probe_1[Ohm]"][1106:1271], p0 = [0.03, 1, 9.5, 0.017])
perr_2 = np.sqrt(np.diag(params_cov_2)) #Error of the fit parameters

print("2 Ampere", params_2[2], "+-", perr_2[2])

#3 Ampere
params_3, params_cov_3 = optimize.curve_fit(Logis, Nb_Si_MF["T[K]"][1428:1573], \
    Nb_Si_MF["R_Probe_1[Ohm]"][1428:1573], p0 = [0.03, 1, 9.5, 0.017])
perr_3 = np.sqrt(np.diag(params_cov_3)) #Error of the fit parameters

print("3 Ampere", params_3[2], "+-", perr_3[2])

#4 Ampere
params_4, params_cov_4 = optimize.curve_fit(Logis, Nb_Si_MF["T[K]"][1705:1870], \
    Nb_Si_MF["R_Probe_1[Ohm]"][1705:1870], p0 = [0.03, 1, 9.5, 0.017])
perr_4 = np.sqrt(np.diag(params_cov_4)) #Error of the fit parameters

print("4 Ampere", params_4[2], "+-", perr_4[2])

#5 Ampere
params_5, params_cov_5 = optimize.curve_fit(Logis, Nb_Si_MF["T[K]"][2061:2261], \
    Nb_Si_MF["R_Probe_1[Ohm]"][2061:2261], p0 = [0.03, 1, 9.5, 0.017])
perr_5 = np.sqrt(np.diag(params_cov_5)) #Error of the fit parameters

print("5 Ampere", params_5[2], "+-", perr_5[2])

#6 Ampere
params_6, params_cov_6 = optimize.curve_fit(Logis, Nb_Si_MF["T[K]"][2481:2761],\
    Nb_Si_MF["R_Probe_1[Ohm]"][2481:2761], p0 = [0.03, 1, 9.5, 0.017])
perr_6 = np.sqrt(np.diag(params_cov_6)) #Error of the fit parameters

print("6 Ampere", params_6[2], "+-", perr_6[2])

#7 Ampere
params_7, params_cov_7 = optimize.curve_fit(Logis, Nb_Si_MF["T[K]"][2906:2995], \
    Nb_Si_MF["R_Probe_1[Ohm]"][2906:2995], p0 = [0.03, 1, 9.5, 0.017])
perr_7 = np.sqrt(np.diag(params_cov_7)) #Error of the fit parameters

print("7 Ampere", params_7[2], "+-", perr_7[2])



#Plot
plt.rc('xtick', labelsize=11)
plt.rc('ytick', labelsize=11)
plt.gcf().subplots_adjust(bottom=0.1)
plt.gcf().subplots_adjust(left=0.15)
plt.gcf().subplots_adjust(top=0.9)


plt.errorbar(Nb_Si_2["T[K]"], Nb_Si_2["R_Probe_1[Ohm]"], fmt='x')
plt.errorbar(Nb_Si_MF["T[K]"][780:920], Nb_Si_MF["R_Probe_1[Ohm]"][780:920], fmt='x')
plt.errorbar(Nb_Si_MF["T[K]"][1106:1271], Nb_Si_MF["R_Probe_1[Ohm]"][1106:1271], fmt='x')
plt.errorbar(Nb_Si_MF["T[K]"][1428:1573], Nb_Si_MF["R_Probe_1[Ohm]"][1428:1573], fmt='x')
plt.errorbar(Nb_Si_MF["T[K]"][1705:1870], Nb_Si_MF["R_Probe_1[Ohm]"][1705:1870], fmt='x')
plt.errorbar(Nb_Si_MF["T[K]"][2061:2261], Nb_Si_MF["R_Probe_1[Ohm]"][2061:2261], fmt='x')
plt.errorbar(Nb_Si_MF["T[K]"][2481:2761], Nb_Si_MF["R_Probe_1[Ohm]"][2481:2761], fmt='x')
plt.errorbar(Nb_Si_MF["T[K]"][2906:2995], Nb_Si_MF["R_Probe_1[Ohm]"][2906:2995], fmt='x')


plt.errorbar(np.linspace(7,11, 1000), Logis(np.linspace(7,11, 1000), \
    params[0], params[1], params[2], params[3]), label='Fitted function', c= 'black')

plt.errorbar(np.linspace(7,11, 1000), Logis(np.linspace(7,11, 1000), \
    params_1[0], params_1[1], params_1[2], params_1[3]), label='Fitted function', c= 'black')

plt.errorbar(np.linspace(7,11, 1000), Logis(np.linspace(7,11, 1000), \
    params_2[0], params_2[1], params_2[2], params_2[3]), label='Fitted function', c= 'black')

plt.errorbar(np.linspace(7,11, 1000), Logis(np.linspace(7,11, 1000), \
    params_3[0], params_3[1], params_3[2], params_3[3]), label='Fitted function', c= 'black')

plt.errorbar(np.linspace(7,11, 1000), Logis(np.linspace(7,11, 1000), \
    params_4[0], params_4[1], params_4[2], params_4[3]), label='Fitted function', c= 'black')

plt.errorbar(np.linspace(7,11, 1000), Logis(np.linspace(7,11, 1000), \
    params_5[0], params_5[1], params_5[2], params_5[3]), label='Fitted function', c= 'black')

plt.errorbar(np.linspace(7, 11, 1000), Logis(np.linspace(7,11, 1000), \
    params_6[0], params_6[1], params_6[2], params_6[3]), label='Fitted function', c= 'black')

plt.errorbar(np.linspace(7, 11, 1000), Logis(np.linspace(7,11, 1000), \
    params_7[0], params_7[1], params_7[2], params_7[3]), label='Fitted function', c= 'black')

plt.xlim(7, 11)
plt.ylim(0.015, 0.03)

plt.title('Super Conductor with Magnetic Field', fontsize = 20)
plt.xlabel('T [K]', fontsize = 13)
plt.ylabel('$R_{SuperCond}$ [$\Omega$]', fontsize = 13)

plt.legend([r'$T_{c, 0A}$ = ' + str(round(params[2], 2)) + " $\pm$ " + str(round(perr[2], 3)) , \
 r'$T_{c, 1A}$ = ' + str(round(params_1[2], 2)) + " $\pm$ " + str(round(perr_1[2], 3)) , \
 r'$T_{c, 2A}$ = ' + str(round(params_2[2], 2)) + " $\pm$ " + str(round(perr_2[2], 3)) , \
 r'$T_{c, 3A}$ = ' + str(round(params_3[2], 2)) + " $\pm$ " + str(round(perr_3[2], 3)) , \
 r'$T_{c, 4A}$ = ' + str(round(params_4[2], 2)) + " $\pm$ " + str(round(perr_4[2], 3)) , \
 r'$T_{c, 5A}$ = ' + str(round(params_5[2], 2)) + " $\pm$ " + str(round(perr_5[2], 3)) , \
 r'$T_{c, 6A}$ = ' + str(round(params_6[2], 2)) + " $\pm$ " + str(round(perr_6[2], 3)) , \
 r'$T_{c, 7A}$ = ' + str(round(params_7[2], 2)) + " $\pm$ " + str(round(perr_7[2], 3))], fontsize = 12)

plt.savefig('Plots/Nb_and_Si_MF.png', dpi=300)
plt.clf()

#Bio Sarvat
#Defining function
def Bio(I):             #Bio Sarvat for the used coil
    return 0.03233 * I

'''There must be an extra factor included due to the fact, that the probe is not placed in
the centre of the coil and therefore does not experience a homogenous magnetic field.
Centre = (33.5cm - 21.5cm)/2 + 21.5cm = 28.5cm, Probe at: 27.5cm --> 1cm. However, the deviation
is so small, that it does not play a dmoniant role.'''

#Read data
Tc_over_I = pd.DataFrame()
Tc_over_I['Current'] = np.linspace(0, 7, 8) #7 different currents
Tc_over_I['Tc'] = [params[2], params_1[2], params_2[2], params_3[2],\
    params_4[2], params_5[2], params_6[2], params_7[2]] # Corresponding critical temperatures
Tc_over_I["B"] = [Bio(0), Bio(1), Bio(2), Bio(3), Bio(4), Bio(5), Bio(6), Bio(7)] #Corresponding Fields

phi = 2.07e-15          #From script
T_c = params[2]         #Critical temperature without magnetic field

#Fit
def Fit(T, zeta):
    return  phi/ (2*math.pi*(zeta)**2)  * (1 - T/T_c)

par, par_cov = optimize.curve_fit(Fit, Tc_over_I["Tc"], Tc_over_I["B"], p0 = 2e-8)  #Fitting critical temperature
err = np.sqrt(np.diag(par_cov)) #Error via fitting parameter

Bc_0 = phi / (2 * np.pi * par[0]**2)   #Calculating Bc(0)
Bc_02 = Bio(1)/(1-8.8/9.35)         #Bc2 close to the critical temperature should yield a better approximation
print("Bc_2 for closest measure point", Bc_02)
errB = -Fit(0, par[0] + err[0]) + Bc_0   #Calculating Bc(0) via Fit: Upper slope straight with y axis


l = par[0]**2/39e-9
lerr = 2 * par[0]/39e-9 * err

#Plot
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.15)
plt.gcf().subplots_adjust(top=0.85)
plt.gcf().subplots_adjust(right=0.85)

plt.errorbar(Tc_over_I["Tc"], Tc_over_I["B"], c= 'black', fmt = "^", markersize='16')
plt.errorbar(np.linspace(7.42, 9.5), Fit(np.linspace(7.5, 9.5), par[0] ))
plt.errorbar(np.linspace(7.42, 9.5), Fit(np.linspace(7.5, 9.5), (par[0] - err[0])))
plt.errorbar(np.linspace(7.42, 9.5), Fit(np.linspace(7.5, 9.5), (par[0] + err[0])) )
plt.xlim(7, 10)
plt.ylim(-0.02, 0.23)

plt.title('Critical magnetic field\n', fontsize = 22)
plt.ylabel('$B_c[T]$', fontsize = 16)
plt.xlabel('$T[K]$', fontsize = 16)

plt.text( 8.0, 0.2, "$\u03BE_{GL}^0$ = (" + str(round(par[0], 10)) + \
    " $\pm$ " + str(round(err[0], 10)) +") m", fontsize = 15 )  #Zeta GL label

plt.text( 8.2, 0.17, "$B_{c}(0)$ = (" + str(round(Bc_0, 2)) + \
    " $\pm$ " + str(round(errB, 2)) +") T", fontsize = 15 )   #B_c2(0) label

plt.text( 8.2, 0.14, "$l$ = (" + str(round(l, 10)) + \
    " $\pm$ " + str(round(lerr[0], 10)) +") m", fontsize = 15 )   #B_c2(0) label




plt.savefig('Plots/Tc_over_I.png', dpi=300)
plt.clf()
