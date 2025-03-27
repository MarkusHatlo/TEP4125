import matplotlib.pyplot as plt
import numpy as np
mole_fraction_table = [0,0.2,0.4,0.5,0.6,0.85,1.0]
molar_Volum_table = [0.09,0.18,0.21,0.22,0.23,0.25,0.255]
mole_frac = 0.6

molar_volum_gradient = np.gradient(molar_Volum_table, mole_fraction_table)

index = mole_fraction_table.index(0.6)

slope = molar_volum_gradient[index]

molal_ch4 = slope * (1 - mole_frac) + molar_Volum_table[index]
molal_CO2 = slope * (0 - mole_frac) + molar_Volum_table[index]



print(f"Slope at mole_fraction 0.6: {slope}")
print(f'The partial molal volume of Ch4 {molal_ch4}')
print(f'The partial molal volume of CO2 {molal_CO2}')


plt.figure(figsize=(10, 8))
plt.plot(mole_fraction_table, molar_Volum_table, label='Molar Volume Curve')

x_tangent = np.linspace(0, 1, 100)
y_tangent = slope * (x_tangent - mole_frac) + molar_Volum_table[index]
plt.plot(x_tangent, y_tangent, '--', label='Tangent at mole_fraction 0.6')

plt.xlabel('Mole Fraction')
plt.ylabel('Molar Volume')
plt.legend()
plt.grid()
plt.show()