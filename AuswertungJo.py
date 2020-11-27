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


#Plot
plt.errorbar(Nb_Si_2["T[K]"], Nb_Si_2["R_Probe_1[Ohm]"], \
    xerr = 0, yerr = 0, c='blue', ecolor='red', fmt='x')
plt.errorbar(np.linspace(7,11, 1000), Logis(np.linspace(7,11, 1000), \
    params[0], params[1], params[2], params[3]), label='Fitted function', c= 'black')

plt.xlim(8, 10.5)
plt.ylim(0.015, 0.0325)

plt.title('Super Conductor without Magnetic Field', fontsize = 20)
plt.xlabel('T [K]', fontsize = 14)
plt.ylabel('$R_{SuperCond}$ [$\Omega$]', fontsize = 14)

plt.text( 8.25, 0.027, "$\Theta_0$ = " + str(round(params[2], 2)) + \
    " $\pm$ " + str(round(perr[2], 3)), fontsize = 15 )
plt.legend(['Measured Resistivity',  "Fitted Plot"], fontsize = 12)

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
plt.xlabel('T [K]')
plt.ylabel('$R_{SuperCond}$ [$\Omega$]')

plt.legend([r'$\Theta_{0, 0A}$ = ' + str(round(params[2], 2)) + " $\pm$ " + str(round(perr[2], 3)) , \
 r'$\Theta_{0, 1A}$ = ' + str(round(params_1[2], 2)) + " $\pm$ " + str(round(perr_1[2], 3)) , \
 r'$\Theta_{0, 2A}$ = ' + str(round(params_2[2], 2)) + " $\pm$ " + str(round(perr_2[2], 3)) , \
 r'$\Theta_{0, 3A}$ = ' + str(round(params_3[2], 2)) + " $\pm$ " + str(round(perr_3[2], 3)) , \
 r'$\Theta_{0, 4A}$ = ' + str(round(params_4[2], 2)) + " $\pm$ " + str(round(perr_4[2], 3)) , \
 r'$\Theta_{0, 5A}$ = ' + str(round(params_5[2], 2)) + " $\pm$ " + str(round(perr_5[2], 3)) , \
 r'$\Theta_{0, 6A}$ = ' + str(round(params_6[2], 2)) + " $\pm$ " + str(round(perr_6[2], 3)) , \
 r'$\Theta_{0, 7A}$ = ' + str(round(params_7[2], 2)) + " $\pm$ " + str(round(perr_7[2], 3))])




plt.savefig('Plots/Nb_and_Si_MF.png', dpi=300)
plt.clf()
