import math
import numpy as np
from firstpart import An, AtAn

#Gas constants
R = 287
ga = 1.4
gg = 1.33
cpa = 1000
cpg = 1150
LHV = 43e6

#Turbine Inlet Temperature
Tt4 = 1150

#Turbine to nozzle inlet area ratio from part 1


#Cruise conditions (AMBIENT TEMP AND PRESSURE)
T_inf = 218.808
P_inf = 23842.3

#Efficiencies
etaNozz = 1
etaComp = 1
etaTurb = 1
etaComb = 1
etaMech = 1
#Pressure ratios
prInlet = 1
prComb = 1

#Exhaust gas mass flow rate
# mg = ma + mf

#Mach number
M = 0.78
v_inf = M*math.sqrt(ga*R*T_inf)

#Stage definition
#---------------------------------
#0 - Atmosphere
#1 - before Inlet
#2 - after Inlet
#3 - after Compressor
#4 - after Combustion
#5 - after Turbine
#6 - after Nozzle

#Inlet
Tt2 = T_inf
Pt2 = P_inf * prInlet

Tt4Tt2 = Tt4/Tt2

print('T04 over T02 is ',Tt4Tt2)

kH = 1 - (AtAn)**((2*ga - 2)/(ga + 1))

OPR = (1 + (cpg/cpa)*kH*Tt4Tt2)**(ga / (ga-1))
print('OPR is', OPR)

mrc = 1.281 * AtAn * OPR * np.sqrt(cpa / cpg / Tt4Tt2)

print('mrc is', mrc)
print('An', An)

#Compressor
Tt3 = Tt2 * (1 + (OPR ** ((ga - 1) / ga) - 1) / etaComp)
Pt3 = Pt2 * OPR
Wcomp = cpa * (Tt3 - Tt2)

#Combustion
#mf = (ma * cpg * (Tt4 - Tt3))/(etaComb * LHV)
Pt4 = Pt3 * prComb

#mg = ma + mf

#Turbine

Tt5 = Tt4 - (Wcomp / cpg)
Pt5 = Pt4 * ((1 - (1 - Tt5 / Tt4) / etaTurb) ** (gg / (gg - 1)))


#print('The pt5_og is {} and the pt5 is {}'.format(Pt5_og, Pt5))

#Nozzle
Pt6 = Pt5
Tt6 = Tt5

critPressRatio = (1 - 1 / etaNozz * ((gg - 1) / (gg +1))) ** (-1 / gg) #Critical pressure ratio

if critPressRatio > Pt6 / P_inf: #Nozzle is unchoked
 print('The nozzle is unchoked')
 v6 = math.sqrt(2 * cpg * etaNozz * Tt6 * (1 - (Pt6 / P_inf) ** ((1 - gg) / gg)))
 Ts6 = Tt6 - (1 / (2 * cpg) * v6 ** 2)
 Ps6 = P_inf
else: #Nozzle is choked
 print('The nozzle is choked')
 Ts6 = Tt6 * (2 / (gg + 1))
 Ps6 = Pt6 / critPressRatio
 v6 = math.sqrt(gg * R * Ts6)

rho6 = Ps6 / (R * Ts6) #Static density at nozzle exit
#A6 = mg / rho6 * v6 #Nozzle exit area

mg = An * rho6 * v6
print('mg is', mg)

vEff6 = v6 + (An * (Ps6 - P_inf)) / mg #Effective velocity at nozzle exit
#Fnew = mg * (vEff6 - v6) #New trust
#Fnew = mg * (vEff6 - v_inf) #New gross trust
Fnew = mg * v6 + An * (Ps6 - P_inf)
Fnew2 = mg * (v6 - v_inf) + An * (Ps6 - P_inf)


print('The new gross thrust is: ', Fnew, 'N')
print('The new net thrust is: ', Fnew2, 'N')

#calculating At/An

#AtAn = (Pt5/Pt4)**((gg+1)/(2*gg))
#print(AtAn)

