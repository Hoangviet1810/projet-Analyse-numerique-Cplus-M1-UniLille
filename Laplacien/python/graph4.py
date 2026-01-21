import numpy as np
import matplotlib.pyplot as plt

T = np.loadtxt("T_final.txt")

plt.figure(figsize=(6,5))
plt.contourf(T, levels=20, cmap="hot")
plt.colorbar(label="Température (°C)")
plt.contour(T, levels=10, colors="black", linewidths=0.5)

plt.title("Lignes de niveau de la température finale T")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.show()
