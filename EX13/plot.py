# quick python one-liner on login node
import numpy as np, matplotlib.pyplot as plt


k,P,_ = np.loadtxt('output_EX6.txt',unpack=True)
plt.loglog(k,P,'.'); 
plt.xlabel('k'); 
plt.ylabel('P(k)'); 
plt.savefig("log-log.png")

