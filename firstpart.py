import math

#Gross thrust
F = 15167
#Air mass flow rate
ma = 23.81
#Fuel mass flow rate
mf = 0.4267
#OPR
OPR = 5.5
#Sea level conditions
T0 = 288
P0 = 1e5
#Gas constants
R = 287
ga = 1.4
gg = 1.33
cpa = 1000
cpg = 1150
LHV = 43e6
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
mg = ma + mf

#Mach number
M = 0
v_inf = M * math.sqrt(ga*R*T0)

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



Tt2 = T0
Pt2 = P0 * prInlet

#Compressor
Tt3 = Tt2 * (1 + (OPR ** ((ga - 1) / ga) - 1) / etaComp)
Pt3 = Pt2 * OPR
Wcomp = cpa * (Tt3 - Tt2)

#Combustion
Tt4 = Tt3 + (mf * LHV * etaComb) / (mg * cpg)
Pt4 = Pt3 * prComb

#Turbine
Tt5 = Tt4 - (Wcomp / cpg)
Pt5 = Pt4 * ((1 - (1 - Tt5 / Tt4) / etaTurb) ** (gg / (gg - 1)))

#Nozzle
Pt6 = Pt5
Tt6 = Tt5

critPressRatio = (1 - 1 / etaNozz * ((gg - 1) / (gg +1))) ** (-1 / gg) #Critical pressure ratio


if critPressRatio > Pt6 / P0: #Nozzle is unchoked
 print('The nozzle is unchoked')
 v6 = math.sqrt(2 * cpg * etaNozz * Tt6 * (1 - (Pt6 / P0) ** ((1 - gg) / gg)))
 Ts6 = Tt6 - (1 / (2 * cpg) * v6 ** 2)
 Ps6 = P0
else: #Nozzle is choked
 print('The nozzle is choked')
 Ts6 = Tt6 * (2 / (gg + 1))
 Ps6 = Pt6 / critPressRatio
 v6 = math.sqrt(gg * R * Ts6)


rho6 = Ps6 / (R * Ts6) #Static density at nozzle exit
A6 = mg / (rho6 * v6) #Nozzle exit area
An = A6

vEff6 = v6 + (A6 * (Ps6 - P0)) / mg #Effective velocity at nozzle exit
#Fnew = mg * (vEff6 - v6) #New trust
#Fnew = mg * (vEff6 - v_inf) #New gross trust
Fnew = mg * v6 + A6 * (Ps6 - P0)

#calculating At/An

AtAn = (Pt5/Pt4)**((gg+1)/(2*gg))


if __name__ == '__main__':
    print('Pressure ratio is {:.2f}'.format(Pt6/P0))
    print('Critical pressure ratio is {:.2f}'.format(critPressRatio))

    print('The exit nozzle area is', A6)
    print('The new gross thrust is: ', Fnew, 'N')
    print('Turbine/Nozzle area ratio is {:.5f}'.format(AtAn))