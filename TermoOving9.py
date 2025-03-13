import sympy as sp

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

# Task 2
print('\nTask 2')