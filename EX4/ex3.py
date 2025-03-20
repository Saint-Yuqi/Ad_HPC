import numpy as np
old_data = np.fromfile("density_old.dat", dtype=np.float32)
new_data = np.fromfile("density.dat", dtype=np.float32)

# Compare element-wise difference
diff = np.abs(old_data - new_data)
print("Max difference:", diff.max())
print("Mean difference:", diff.mean())

