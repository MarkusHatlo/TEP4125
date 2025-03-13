import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Task 1
print('Task 1')
v_molar = sp.Symbol('v_molar')
p = 20e6
T = 673
a = 5.531e5
b = 0.0305
R = 8.314e3

#van der Waals
eq_vdw = sp.Eq(p, R*T/(v_molar-b) - a/v_molar**2)
v_molar = sp.solve(eq_vdw, v_molar)

v_mol = v_molar[0]
v_spesific = v_mol/18
print('v_spesific using van der Waals:', v_spesific)

#redlich-kwong
a_rk = 142.59e5
b_rk = 0.02111
v_molar_rk = sp.Symbol('v_molar_rk')

eq_rk = sp.Eq(p, R*T/(v_molar_rk-b_rk) - a_rk/(v_molar_rk*(v_molar_rk+b_rk)*T**0.5))
v_molar_rk = sp.solve(eq_rk, v_molar_rk)

v_mol_rk = v_molar_rk[0]
v_spesific_rk = v_mol_rk/18
print('v_spesific using redlich-kwong:', v_spesific_rk)

#Task 3
print('\nTask 3')

R_spes = 8.314e3/44.01
a_Walls_CO2 = 3.647e5/44.01**2
b_Walls_CO2 = 0.0428/44.01**2
a_Redlich_CO2 = 64.43e5/44.01**2
b_Redlich_CO2 = 0.02963/44.01**2
T3 = 240

def calcIdealGas(v):
    return (R_spes*T3/v)/1e5

def calcVanDerWaals(v):
    return ((R_spes*T3/(v-b_Walls_CO2)) - (a_Walls_CO2/v**2))/1e5

def calcRedlichKwong(v):
    return ((R_spes*T3/(v-b_Redlich_CO2)) - a_Redlich_CO2/(v*(v+b_Redlich_CO2)*T3**0.5))/1e5



v_plot = np.linspace(1e-4,0.05,100)
p_Ideal_plot = calcIdealGas(v_plot)
p_VDW_plot = calcVanDerWaals(v_plot)
p_RK_plot = calcRedlichKwong(v_plot)

plt.plot(v_plot, p_Ideal_plot)
plt.plot(v_plot, p_VDW_plot)
plt.plot(v_plot, p_RK_plot)
plt.legend(['Ideal gas', 'Van der Waals', 'Redlich-Kwong'])
plt.ylim(0, 50)
#plt.xlim(0.01, 0.05)
plt.xlabel('Volume [m^3/mol]')
plt.ylabel('Pressure [Bar]')
plt.title('Ideal gas')
plt.show()

