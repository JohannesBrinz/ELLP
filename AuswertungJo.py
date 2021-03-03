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
from scipy import constants
import math
import matplotlib.font_manager as fm

# Cu_Si_1
#Datenimport
Cu_Si_1 = pd.read_csv("Cu_and_Si_1.dat", sep = "\t", header = 4, \
    names = ["Zeit [s]","T [K]", "R_Cu [Ohm]", "R_Thermom [Ohm]", "R_Si [Ohm]"])
Cu_Si_2 = pd.read_csv("Cu_and_Si_2.dat", sep = "\t", header = 4, \
    names = ["Zeit [s]","T [K]", "R_Cu [Ohm]", "R_Thermom [Ohm]", "R_Si [Ohm]"])
Cu_Si_3 = pd.read_csv("Cu_and_Si_3.dat", sep = "\t", header = 4, \
    names = ["Zeit [s]","T [K]", "R_Cu [Ohm]", "R_Thermom [Ohm]", "R_Si [Ohm]"])
#Sort
Cu_Si_1 = Cu_Si_1.sort_values(by = "Zeit [s]")

#Shift Temperature
Cu_Si_3_shifted = []
for i in range(3572) :
     Cu_Si_3_shifted =  Cu_Si_3["T [K]"] - 20

Cu_Si_2_shifted = []
for i in range(90) :
     Cu_Si_2_shifted =  Cu_Si_2["T [K]"] - 20
#COPPER
#Shift Resistance
R_T_1 = []
for i in range(1161) :
     R_T_1 =  Cu_Si_1["R_Cu [Ohm]"] - 0.0169

R_T_2 = []
for i in range(90) :
     R_T_2 =  Cu_Si_2["R_Cu [Ohm]"] - 0.0169

R_T_3 = []
for i in range(3572) :
     R_T_3 =  Cu_Si_3["R_Cu [Ohm]"] - 0.0169


#Calculate RRR

RRR = Cu_Si_3["R_Cu [Ohm]"][3555]/Cu_Si_1["R_Cu [Ohm]"][4]
print("RRR=",RRR)
#Fit
def Fit(T, R_T, T_Debye):
    return (1.17* R_T *T / T_Debye - 0.17*R_T)

param, params_cov = optimize.curve_fit(Fit, Cu_Si_3_shifted, R_T_3 )
print("R_T = ",param[0] )
print("T_Debye = ", param[1])
#Plot Resistance

plt.errorbar(Cu_Si_1["T [K]"][0:788]/param[1], R_T_1[0:788]/param[0], c='red', ecolor='red', fmt='.')

plt.errorbar(Cu_Si_2_shifted[50:93]/param[1], R_T_2[50:93]/param[0], c='blue', ecolor='red', fmt='.')

plt.errorbar(Cu_Si_3_shifted/param[1], R_T_3/param[0], c='orange', ecolor='red', fmt='.')

plt.errorbar(np.linspace(5,300, 10)/param[1], Fit(np.linspace(5,300, 10), param[0], param[1])/param[0], c= 'black')
plt.title('Copper Resistivity', fontsize = 16)
plt.xlabel(R'$\frac{T}{\Theta_D}$', fontsize = 12)
plt.ylabel(R'$\frac{R}{R(\Theta_D)}$', fontsize = 12)
plt.text(0.4, 0.8, "$\Theta_D=$" + str(round(param[1], 1)) + "K", fontsize = 12)
plt.legend(['experimental data termometer 1', 'experimental data termometer 2', 'experimental data termometer 3', R"Grüneisen fit: $  R_T (T) = 1.17 R_T (\Theta_D )\frac{T}{\Theta_D} - 0.17R_T (\Theta_D )$"], fontsize = 10)
plt.savefig('Plots/Cu_and_Si_1.png', dpi=300)
plt.clf()

#Universal Reduced Resistivity Plot
T_red = Cu_Si_1["T [K]"][20:788]/param[1]
R_red = (R_T_1[20:788]-R_T_1[20])/param[0]

print("\nhier:  ", R_T_1[20], R_T_1[788])

#Description fit from Gerthesen
Cu = pd.DataFrame()
Cu["T"] = [0.266, 0.247, 0.082, 0.062]
Cu["R"] = [0.132, 0.113, 0.003, 0.002]

Au = pd.DataFrame()
Au["T"] = [0.328, 0.120, 0.105]
Au["R"] = [0.207, 0.017, 0.009]

Na = pd.DataFrame()
Na["T"] = [0.388, 0.282]
Na["R"] = [0.267, 0.152]

def Grueneisen(T_red, alpha, beta):
    return(  beta * ( 1 / ( (1/T_red)**2 + alpha * (1/T_red) ) )  )
param1, param1_cov = optimize.curve_fit(Grueneisen, T_red, R_red)

print("beta=" , param1[1])
print("alpha= ", param1[0])
plt.errorbar(T_red, R_red)

plt.errorbar(np.linspace(0.025 ,0.27, 1000), Grueneisen(np.linspace(0.025 , 0.27, 1000), param1[0], param1[1]), label='Fitted Function', c= 'black')
plt.errorbar(Cu["T"], Cu["R"], fmt = "s")
plt.errorbar(Au["T"], Au["R"], fmt = "o")
plt.errorbar(Na["T"], Na["R"], fmt = "D")

plt.title('Universal Plot according to Grüneisen Theory', fontsize = 16)
plt.xlabel('Reduced Temperature [-]', fontsize = 12)
plt.ylabel('Reduced Resistance [-]', fontsize = 12)

plt.legend(["Experimental data Cu", "Grüneisen fit", "Literature value Cu", "Literature value Au", "Literature value Na"], fontsize = 12)

plt.savefig('Plots/univ.png', dpi=300)
plt.clf()


#SILICON
#Plot Resistance

plt.errorbar(Cu_Si_1["T [K]"][400:788], Cu_Si_1["R_Si [Ohm]"][400:788] )

plt.errorbar(Cu_Si_2_shifted[50:93], Cu_Si_2["R_Si [Ohm]"][50:93])

plt.errorbar(Cu_Si_3_shifted, Cu_Si_3["R_Si [Ohm]"])

plt.yscale("log")
plt.xlabel("T [K]")
plt.ylabel('R [Ohm]')
plt.title("Resistivity semiconducting doped silicon")

plt.savefig('Plots/Si-Resistance.png', dpi=300)
plt.clf()

#Fit determination E_g

def E_g(T, C, E_g):
    k = constants.Boltzmann
    return C * np.exp(-(E_g)/(2*k*T))

l = 5e-3
A = 9e-3 * 0.5e-3

param2, param2_cov = optimize.curve_fit(E_g, Cu_Si_1["T [K]"][450:500], (l/A)/Cu_Si_1["R_Si [Ohm]"][450:500], p0 = [1e-6, 1e-20])
err2 = np.sqrt(np.diag(param2_cov))

print("Die Parameter sind: ", param2)
print("Und die Fehler: ", err2)
print("Und in meV: ", err2[1]/constants.e)

plt.errorbar(1/Cu_Si_1["T [K]"][450:788], (l/A)/Cu_Si_1["R_Si [Ohm]"][450:788] )

plt.errorbar(1/Cu_Si_2_shifted[50:93], (l/A)/Cu_Si_2["R_Si [Ohm]"][50:93])

plt.errorbar(1/Cu_Si_3_shifted, (l/A)/Cu_Si_3["R_Si [Ohm]"])

plt.errorbar(1/np.linspace(25, 75, 1000), E_g(np.linspace(25, 75, 1000), param2[0],  param2[1]), c="black")

plt.yscale("log")

plt.text( 0.027, 1, "$E_d$ = " + str(round(param2[1]*1000/constants.e )) + "meV", fontsize = 12)
plt.xlabel(r"$\frac{1}{T} [K^{-1}]$", fontsize = 12)
plt.ylabel('$\sigma [m/\Omega]$', fontsize = 12)
plt.title("Conductivity semiconducting doped silicon", fontsize = 16)
plt.legend(['experimental data termometer 1', 'experimental data termometer 2','experimental data termometer 3', "fit"], fontsize = 12)
plt.savefig('Plots/Si-ResistanceII.png', dpi=300)
plt.clf()

#Plot Resistance high T
plt.errorbar(Cu_Si_1["T [K]"][600:788], Cu_Si_1["R_Si [Ohm]"][600:788] )
plt.errorbar(Cu_Si_3_shifted, Cu_Si_3["R_Si [Ohm]"])
plt.xlabel('T [K]')
plt.ylabel('R [Ohm]')
plt.title("R(T)")

plt.savefig('Plots/Si-Resistance-cut.png', dpi=300)
plt.clf()


#Natural log of Conductivity in dependence of inverse T
k_B = constants.k
eV = constants.e
l = 5e-3
A = 9e-3 * 0.5e-3
sigma1 = l/(A*Cu_Si_1["R_Si [Ohm]"][400:788])
lnrsigma1 = np.log(sigma1)
sigma2 = l/(A*Cu_Si_2["R_Si [Ohm]"][60:93])
lnrsigma2 = np.log(sigma2)
sigma3 = l/(A*Cu_Si_3["R_Si [Ohm]"])
lnrsigma3 = np.log(sigma3)

T_inv_1 = 1/ Cu_Si_1["T [K]"][400:788]
T_inv_2 = 1/ Cu_Si_2_shifted[60:93]
T_inv_3 = 1/ Cu_Si_3_shifted

plt.errorbar(T_inv_1, lnrsigma1)

plt.errorbar(T_inv_2, lnrsigma2)

plt.errorbar(T_inv_3, lnrsigma3)

plt.xlabel('1/T [K⁻¹]')
plt.ylabel('log(σ) [$ln((m\Omega)^{-1})$]')
plt.title("Natural log of Conductivity in dependence of inverse Temperature")

plt.savefig('Plots/Si-lnsigma.png', dpi=300)
plt.clf()


#Energy Fit
plt.errorbar(T_inv_3[69:416], lnrsigma3[69:416])
def E_fit(T_inv, E_d, q):
    return(-E_d*T_inv/(2*k_B) + q)
p, p_cov = optimize.curve_fit(E_fit, T_inv_3[69:416], lnrsigma3[69:416])
plt.errorbar(np.linspace(0.0082 ,0.01239, 1000), E_fit(np.linspace(0.0082 , 0.01239, 1000), p[0], p[1]), label='Fitted Function', c= 'black')
plt.xlabel('1/T [K⁻¹]')
plt.ylabel('log(σ) [$ln((m\Omega)^{-1})$]')
plt.title("Natural log of Conductivity in dependence of inverse Temperature")
plt.savefig('Plots/Si-extrinsic-range.png', dpi=300)

E_d_err = np.sqrt(np.diag(p_cov))
print("E_d = ", p[0]/eV, "eV")
print("Delta E_d = ", E_d_err[0]/eV ,"eV")



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

#Bio Sarvat
#Defining function
def Bio(I):             #Bio Sarvat for the used coil
    return 0.03233 * I

'''There must be an extra factor included due to the fact, that the probe is not placed in
the centre of the coil and therefore does not experience a homogenous magnetic field.
Centre = (33.5cm - 21.5cm)/2 + 21.5cm = 28.5cm, Probe at: 27.5cm --> 1cm. However, the deviation
is so small, that it does not play a dmoniant role.'''

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
#plt.grid(True)
plt.legend([r'$T_{c, 0mT}$ = ' + str(round(params[2], 2)), \
 r'$T_{c, ' + str(round(Bio(1)*1000)) + "mT}$="  + str(round(params_1[2], 2)),  \
 r'$T_{c, ' + str(round(Bio(2)*1000)) + "mT}$ = " + str(round(params_2[2], 2)), \
 r'$T_{c, ' + str(round(Bio(3)*1000)) + "mT}$ = " + str(round(params_3[2], 2)), \
 r'$T_{c, ' + str(round(Bio(4)*1000)) + "mT}$ = " + str(round(params_4[2], 2)), \
 r'$T_{c, ' + str(round(Bio(5)*1000)) + "mT}$ = " + str(round(params_5[2], 2)), \
 r'$T_{c, ' + str(round(Bio(6)*1000)) + "mT}$ = " + str(round(params_6[2], 2)), \
 r'$T_{c, ' + str(round(Bio(7)*1000)) + "mT}$ = " + str(round(params_7[2], 2))], fontsize = 12)

plt.savefig('Plots/Nb_and_Si_MF.png', dpi=300)
plt.clf()


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

plt.text( 8.2, 0.14, "$l$ = (" + str(round(l, 10)) + \
    " $\pm$ " + str(round(lerr[0], 10)) +") m", fontsize = 15 )   #B_c2(0) label




plt.savefig('Plots/Tc_over_I.png', dpi=300)
plt.clf()


#RRR Nb

RRR = Nb_Si_2["R_Probe_1[Ohm]"][10]/Nb_Si_2["R_Probe_1[Ohm]"][2700]

print("RRR Nb:   ", RRR)


#Debye Nb
def Fit(T, R_T, T_Debye):
    return (1.17* R_T *T / T_Debye - 0.17*R_T)

param, params_cov = optimize.curve_fit(Fit, Nb_Si_2["T[K]"][10:2000], Nb_Si_2["R_Probe_1[Ohm]"][10:2000]-Nb_Si_2["R_Probe_1[Ohm]"][2000] )
print("R_T = ",param[0] )
print("T_Debye = ", param[1])
#Plot Resistance

plt.errorbar(Nb_Si_2["T[K]"][10:2000]/param[1], (Nb_Si_2["R_Probe_1[Ohm]"][10:2000]-Nb_Si_2["R_Probe_1[Ohm]"][2000])/param[0], c='red', ecolor='red', fmt='.')


plt.errorbar(np.linspace(135,275, 10)/param[1], Fit(np.linspace(135,275, 10), param[0], param[1])/param[0], c= 'black')
plt.title('Niobium Resistivity', fontsize = 16)
plt.xlabel(R'$\frac{T}{\Theta_D}$', fontsize = 12)
plt.ylabel(R'$\frac{R}{R(\Theta_D)}$', fontsize = 12)
plt.text(0.17, 0.1, "$\Theta_D=$" + str(round(param[1], 1)) + "K", fontsize = 12)
plt.legend(['experimental data termometer 1', R"Grüneisen fit: $  R_T (T) = 1.17 R_T (\Theta_D )\frac{T}{\Theta_D} - 0.17R_T (\Theta_D )$"], fontsize = 10)
plt.savefig('Plots/Nb_and_Si_1.png', dpi=300)
plt.clf()
