import matplotlib.pyplot as plt
import numpy as np

# Grid resolutions
resolutions = np.array([100, 200, 500])

# Timing data for each phase:
# Reading times in seconds
reading_time = np.array([4.2208e-05, 8.4416e-05, 4.3542e-05])
# Mass assignment times in seconds
mass_assignment_time = np.array([0.026328, 0.0334267, 0.0676093])
# Projection times in seconds
projection_time = np.array([0.000464, 0.00433371, 0.0579554])

# Create a figure with three subplots for comparison
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

# Linear scale plot
axs[0].plot(resolutions, reading_time, label='Reading Time', marker='o')
axs[0].plot(resolutions, mass_assignment_time, label='Mass Assignment', marker='o')
axs[0].plot(resolutions, projection_time, label='Projection Time', marker='o')
axs[0].set_title("Linear Scale")
axs[0].set_xlabel("Grid Resolution (nGrid)")
axs[0].set_ylabel("Time (s)")
axs[0].legend()
axs[0].grid(True)

# Log-log plot
axs[1].loglog(resolutions, reading_time, label='Reading Time', marker='o')
axs[1].loglog(resolutions, mass_assignment_time, label='Mass Assignment', marker='o')
axs[1].loglog(resolutions, projection_time, label='Projection Time', marker='o')
axs[1].set_title("Log-Log")
axs[1].set_xlabel("Grid Resolution (nGrid, log scale)")
axs[1].set_ylabel("Time (s) [Log Scale]")
axs[1].legend()
axs[1].grid(True, which='both')

# Semilog plot (y-axis logarithmic)
axs[2].semilogy(resolutions, reading_time, label='Reading Time', marker='o')
axs[2].semilogy(resolutions, mass_assignment_time, label='Mass Assignment', marker='o')
axs[2].semilogy(resolutions, projection_time, label='Projection Time', marker='o')
axs[2].set_title("Log-Linear")
axs[2].set_xlabel("Grid Resolution (nGrid)")
axs[2].set_ylabel("Time (s) [Log Scale]")
axs[2].legend()
axs[2].grid(True, which='both')



plt.suptitle("Timing vs Grid Resolution for Each Phase")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

plt.savefig("result_time_resolution.png")
plt.show()
