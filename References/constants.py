### Fundamental constants
k = 8.617*10**(-5) # Boltzmann constant
h = 4.135667696*10**(-15)
e = 1.602*10**(-19)

### State functions
T = 298.15

### ZPE

# Adsorbed
ZPE_H2O = 0.67
ZPE_O = 0.07
ZPE_OH = 0.38
ZPE_OOH = 0.43
ZPE_H = 0.17

# Gas
ZPE_H2O_GAS = 0.57
ZPE_H2_GAS = 0.26

# Combined
ZPE_O_H2 = ZPE_O + ZPE_H2_GAS
ZPE_OH_H = ZPE_OH + 0.5*ZPE_H2_GAS
ZPE_OOH_3H = ZPE_OOH + (3/2)*ZPE_H2_GAS
ZPE_H_H = ZPE_H + 0.5*ZPE_H2_GAS

### TS

# Gas
TS_H2O_GAS = 0.59
TS_H2_GAS = 0.41

# Adsorbed
TS_O = 0
TS_OH = 0
TS_OOH = 0
TS_H = 0
TS_H2O = 0

# Combined
TS_O_H2 = TS_O + TS_H2_GAS
TS_OH_H = TS_OH + 0.5*TS_H2_GAS
TS_OOH_3H = TS_OOH + (3/2)*TS_H2_GAS 
TS_H_H = TS_H + 0.5*TS_H2_GAS

### Solvation energies
SV_OH = -0.4
#SV_OH = 0
SV_OOH = -0.4
SV_H = -0.2
SV_O = 0
SV_H2O = -0.4 # From heat of vaporization
#SV_H2O = 0
SV_OO = 0 # unknown

### Enthalpies
E_H2O_GAS = -14.175
E_H2_GAS = -6.992

### Thermal energies

# Gas
TH_H2_GAS = 0.08
TH_H2O_GAS = 0.12

# Adsorbed
TH_O = 0.04
TH_OH = 0.08
TH_OOH = 0.12
TH_H2O = 0.12


### Construct free energies

# Gas
G_H2O_GAS = E_H2O_GAS + ZPE_H2O_GAS + TH_H2O_GAS - TS_H2O_GAS
G_H2_GAS = E_H2_GAS + ZPE_H2_GAS + TH_H2_GAS - TS_H2_GAS

# Adsorbed

G_H2O = ZPE_H2O + TH_H2O - TS_H2O
G_O_H2 = E_H2_GAS + ZPE_O_H2 + TH_O + TH_H2_GAS - TS_O_H2
G_OH_H = 0.5*E_H2_GAS + ZPE_OH_H + TH_OH + 0.5*TH_H2_GAS - TS_OH_H
G_OOH_3H = (3/2)*E_H2_GAS + ZPE_OOH_3H + TH_OOH + (3/2)*TH_H2_GAS - TS_OOH_3H
G_OO_4H = 2*(E_H2_GAS +TH_H2_GAS) # Needs O2 vib
G_H_H = 0.5*(E_H2_GAS+TH_H2_GAS) + ZPE_H_H - TS_H_H





