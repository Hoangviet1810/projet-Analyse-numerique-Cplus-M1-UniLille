import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))

# PCG
pcg_ci = np.loadtxt("residus_PCG_CI.txt")
plt.semilogy(pcg_ci[:,0], pcg_ci[:,1], 'b-o',
             linewidth=2, markersize=4,
             label="PCG : C = I")

pcg_diag = np.loadtxt("residus_PCG_diag.txt")
plt.semilogy(pcg_diag[:,0], pcg_diag[:,1], 'g-s',
             linewidth=2, markersize=4,
             label="PCG : C = diag(A)")

pcg_ic0 = np.loadtxt("residus_PCG_IC0_omega_0.txt")
plt.semilogy(pcg_ic0[:,0], pcg_ic0[:,1], 'r-^',
             linewidth=2, markersize=5,
             label=r"PCG : IC(0), $\omega = 0$")

# (optionnel) IC(0) avec omega = 1/n
try:
    pcg_ic0_w = np.loadtxt("residus_PCG_IC0_omega_1n.txt")
    plt.semilogy(pcg_ic0_w[:,0], pcg_ic0_w[:,1], 'm-d',
                 linewidth=2, markersize=5,
                 label=r"PCG : IC(0), $\omega = 1/n$")
except:
    pass

plt.xlabel("Iteration k")
plt.ylabel(r"$\|r_k\| / \|b\|$")
plt.title("PCG – Cas (1.1) : rho = 1  --> Laplacien standard, N = 20")
plt.legend(fontsize=9)
plt.tight_layout()
plt.show()

