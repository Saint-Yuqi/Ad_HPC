import numpy as np
import matplotlib.pyplot as plt

nGrid = 100
data = np.fromfile("projected_density_PCS_100.dat", dtype=np.float32).reshape((nGrid, nGrid))
plt.imshow(data, cmap="inferno", origin="lower")
plt.colorbar(label="Density")
plt.xlabel("X Grid Index")
plt.ylabel("Y Grid Index")
plt.title("Projected Dark Matter Density (Grid=100)")
plt.savefig("result_100_PCS.png")
plt.close()


# nGrid = 200
# data = np.fromfile("projected_density_200.dat", dtype=np.float32).reshape((nGrid, nGrid))
# plt.imshow(data, cmap="inferno", origin="lower")
# plt.colorbar(label="Density")
# plt.xlabel("X Grid Index")
# plt.ylabel("Y Grid Index")
# plt.title("Projected Dark Matter Density (Grid=200)")
# plt.savefig("result_200.png")
# plt.close()


# nGrid = 500
# data = np.fromfile("projected_density_500.dat", dtype=np.float32).reshape((nGrid, nGrid))
# plt.imshow(data, cmap="inferno", origin="lower")
# plt.colorbar(label="Density")
# plt.xlabel("X Grid Index")
# plt.ylabel("Y Grid Index")
# plt.title("Projected Dark Matter Density (Grid=500)")
# plt.savefig("result_500.png")


