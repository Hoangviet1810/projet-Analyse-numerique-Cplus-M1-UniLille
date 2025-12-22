import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))

# ---------- Steepest descent ----------
sd = np.loadtxt("residus_steepest_descent.txt")
plt.semilogy(sd[:,0], sd[:,1], 'k-', linewidth=2,
             label="Steepest descent")

# ---------- Gradient pas fixe ----------
for j in range(1, 7):
    data = np.loadtxt(f"residus_pas_fixe_j{j}.txt")
    alpha_label = f"Pas fixe j={j}"
    plt.semilogy(data[:,0], data[:,1], '--',
                 label=alpha_label)

plt.xlabel("Iteration k")
plt.ylabel(r"$||r_k|| / ||b||$")
plt.title("Comparaison steepest descent / gradient à pas fixe")
plt.legend()
#plt.grid(True, which="both")

plt.tight_layout()
plt.show()

