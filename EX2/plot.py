import numpy as np
import matplotlib.pyplot as plt

# Load projected density map
nGrid = 100  # Set same grid size
data = np.fromfile("projected_density.dat", dtype=np.float32).reshape((nGrid, nGrid))

# Plot
plt.imshow(data, cmap="inferno", origin="lower")
plt.colorbar(label="Density")
plt.xlabel("X Grid Index")
plt.ylabel("Y Grid Index")
plt.title("Projected Dark Matter Density")
plt.savefig("result.png")
plt.show()
