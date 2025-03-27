import numpy as np
import matplotlib.pyplot as plt

# List your filenames here:
filenames = [
    "power_log_100.dat",
    "power_log_200.dat",
    "power_log_500.dat"
]

# Create a figure with 3 subplots stacked vertically
fig, axes = plt.subplots(nrows=3, figsize=(6, 10))

for i, fname in enumerate(filenames):
    # Load the data from file, ignoring lines that begin with '#'
    data = np.loadtxt(fname)
    # Suppose columns are: [binIndex, k_center, P(k)]
    kcenter = data[:, 1]  # second column
    Pk      = data[:, 2]  # third column

    axes[i].loglog(kcenter, Pk, marker='o', label=fname)
    axes[i].set_xlabel("k")
    axes[i].set_ylabel("P(k)")
    axes[i].set_title(f"Log-binned Power Spectrum: {fname}")
    axes[i].grid(True)
    axes[i].legend()

fig.tight_layout()
plt.savefig("three_subplots.png")
plt.show()  # remove if running on HPC without display
