from firstpart import *
import numpy as np
from secondpart import OPR, mg


''' _____________Step 1_____________'''

#T_inf = 218.808
#P_inf = 23842.3

T0 = 218.808 # [K]
P0 = 23842.3 # [Pa]
Rspecific = 287

# mg = 12.371750426650197
# Omega is selected subject to rotational velocity contraint of 600 m/s
# Omega = 3500

beta_tt = OPR

''' _____________Step 2_____________'''

eta_p = 1# 0.9 # polytropic efficiency


Tt2 = T0

Tt3 = Tt2 * (beta_tt)**((ga-1)/(ga*eta_p))
print('Tt3 is {:.2f} K'.format(Tt3))

w_sp = (Tt3 - Tt2)*cpa
print('Specific Work is {:.2f} J'.format(w_sp))

# Selecting work and flow coefficients and degree of reaction

''' _____________Step 3_____________'''

psi     = 0.25
phi     = 0.3
R       = 0.5

alpha_1 = np.arctan((-R-(psi/2)+1)/(phi))
beta_1  = np.arctan(np.tan(alpha_1) - (1/phi))
beta_2  = np.arctan((psi+(phi*np.tan(alpha_1))-1)/phi)
alpha_2 = np.arctan(np.tan(beta_2) + (1/phi))

print('\n')
print('alpha_1 is {:.2f}°'.format(np.rad2deg(alpha_1)))
print('alpha_2 is {:.2f}°'.format(np.rad2deg(alpha_2)))
print('beta_1 is {:.2f}°'.format(np.rad2deg(beta_1)))
print('beta_2 is {:.2f}°'.format(np.rad2deg(beta_2)))


''' _____________Step 4_____________'''


# Selecting number of stages
stageno = 8

# finding specific work per stage
w_stage = w_sp/stageno


U_stage = np.sqrt(w_stage/psi) * np.ones(stageno + 1)
print('\nU_stage is ',U_stage)


''' _____________Step 5_____________'''

r_stage = np.linspace(0.12,0.12,stageno)
# U_stage = r_stage * Omega
# print(U_stage)


V_1 = phi * U_stage / np.cos(alpha_1)
V_m = V_1 * np.cos(alpha_1)
V_2 = phi * U_stage / np.cos(alpha_2)
W_1 = V_1 * np.cos(alpha_1)/np.cos(beta_1)
W_2 = V_2 * np.cos(alpha_2)/np.cos(beta_2)

print('\n')
print('V1 is {:.2f} m/s'.format(V_1[0]))
print('V2 is {:.2f} m/s'.format(V_2[0]))
print('W1 is {:.2f} m/s'.format(W_1[0]))
print('W2 is {:.2f} m/s'.format(W_2[0]))

''' _____________Step 6_____________'''


# Finding total temperature between each stage
deltaT = psi*(U_stage**2)/cpa

T = T0 * np.ones(len(r_stage)+1)
for i,t in enumerate(T):
    T[i] = T[i] + np.sum(deltaT[0:i])

print('\nTotal temperature between each stage is \n', T)

# Finding total pressure between each stage
p = np.zeros(stageno+1)
p[0] = P0

for i in range(stageno):
    p_new = p[i] * (T[i+1]/T[i])**(eta_p*ga/(ga-1))
    p[i+1] = p_new

print('\nTotal pressure between each stage is \n', p)

rho = p / (Rspecific * T)

print('\nTotal density between each stage is \n', rho)

Area = mg / (rho * V_m)
print('\nArea between each stage is \n', Area)

h = Area / (1 * np.pi * r_stage[0])
print('\nHeight between each stage is \n', h)




