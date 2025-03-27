import numpy as np
import matplotlib.pyplot as plt

# Constants for CO2
T = 240  # Temperature in K
R_spes = 8.314 / 44.01  # Gas constant in kJ/(kg·K)

# Van der Waals constants for CO2
a_vdw = 3.647
b_vdw = 0.0428

# Redlich-Kwong constants for CO2
a_rk = 64.43
b_rk = 0.02963

# Benedict-Webb-Rubin constants for CO2
a_bwr = 0.1386
A_bwr = 2.7737
b_bwr = 0.007210
B_bwr = 0.04991
c_bwr = 1.512e4
C_bwr = 1.404e5
alfa_bwr = 8.47e-5
gamma_bwr = 0.00539

# Define equations of state
def calcIdealGas(v):
    return (R_spes * T * 1000 / v) / 1e5  # Convert kJ to J and Pa to bar

def calcVanDerWaals(v):
    return ((R_spes * T * 1000 / (v - b_vdw)) - (a_vdw / v**2)) / 1e5

def calcBerthelot(v):
    # Berthelot equation using van der Waals constants but with temperature factor
    return ((R_spes * T * 1000 / (v - b_vdw)) - (a_vdw / (T * v**2))) / 1e5

def calcRedlichKwong(v):
    return ((R_spes * T * 1000 / (v - b_rk)) - a_rk / (v * (v + b_rk) * T**0.5)) / 1e5

def calcBenedictWebbRubin(v):
    # Benedict-Webb-Rubin equation
    term1 = R_spes * T * 1000 / v
    term2 = (B_bwr * R_spes * T * 1000 - A_bwr - C_bwr / T**2) / v**2
    term3 = (b_bwr * R_spes * T * 1000 - a_bwr) / v**3
    term4 = a_bwr * alfa_bwr / v**6
    term5 = c_bwr / (v**3 * T**2) * (1 + gamma_bwr / v**2) * np.exp(-gamma_bwr / v**2)
    return (term1 + term2 + term3 + term4 + term5) / 1e5

# Create a range for specific volume (m³/kg)
v_low = np.linspace(0.0006, 0.01, 500)  # More points at low volumes
v_high = np.linspace(0.01, 0.05, 200)   # Fewer points at high volumes
v_plot = np.unique(np.concatenate([v_low, v_high]))

# Calculate pressure for each equation
p_Ideal_plot = calcIdealGas(v_plot)
p_VDW_plot = calcVanDerWaals(v_plot)
p_Berthelot_plot = calcBerthelot(v_plot)
p_RK_plot = calcRedlichKwong(v_plot)
p_BWR_plot = calcBenedictWebbRubin(v_plot)

# Create the plot
plt.figure(figsize=(12, 8))

# Plot each equation
plt.plot(v_plot, p_Ideal_plot, 'b-', label='Ideal Gas Law')
plt.plot(v_plot, p_VDW_plot, 'orange', label='Van der Waals')
plt.plot(v_plot, p_Berthelot_plot, 'g-', label='Berthelot')
plt.plot(v_plot, p_RK_plot, 'r-', label='Redlich-Kwong')
plt.plot(v_plot, p_BWR_plot, 'purple', label='Benedict-Webb-Rubin')

# Add saturation points
p_sat = 12.85  # bar
v_f = 0.000918  # m³/kg (liquid phase)
v_g = 0.02997  # m³/kg (vapor phase)
plt.scatter([v_f, v_g], [p_sat, p_sat], color='red', marker='o', s=100, label='Saturation Points')

# Add horizontal line at saturation pressure
plt.axhline(y=p_sat, color='r', linestyle='--', alpha=0.5)

# Add vertical lines at saturation volumes
plt.axvline(x=v_f, color='r', linestyle='--', alpha=0.3)
plt.axvline(x=v_g, color='r', linestyle='--', alpha=0.3)

# Set plot properties
plt.xlim(0, 0.05)
plt.ylim(0, 50)
plt.grid(True)
plt.xlabel('Specific Volume (m³/kg)', fontsize=14)
plt.ylabel('Pressure (bar)', fontsize=14)
plt.title('Pressure vs. Specific Volume for CO2 at T = 240 K', fontsize=16)
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()
