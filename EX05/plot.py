import numpy as np
import matplotlib.pyplot as plt

# 1) Plot linear-binned data
linear_data = np.loadtxt('power_linear.dat')
k_lin = linear_data[:, 0]
P_lin = linear_data[:, 1]

plt.figure()  # create a new figure
plt.loglog(k_lin, P_lin, label='Linear Binning')
plt.xlabel('k')
plt.ylabel('P(k)')
plt.title('Power Spectrum - Linear Binning')
plt.legend()
plt.savefig('linear_binning.png')  # save to PNG
plt.show()  # or remove this line if you're on an HPC cluster without a display

# 2) Plot variable-binned data
var_data = np.loadtxt('power_variable_3.dat')
k_var = var_data[:, 0]
P_var = var_data[:, 1]

plt.figure()
plt.loglog(k_var, P_var, color='tab:orange', label='Variable Binning')
plt.xlabel('k')
plt.ylabel('P(k)')
plt.title('Power Spectrum - Variable Binning')
plt.legend()
plt.savefig('variable_binning.png')
plt.show()
