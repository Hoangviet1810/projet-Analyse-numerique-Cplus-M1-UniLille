import numpy as np
import matplotlib.pyplot as plt

def plot_erreur(file_path="erreur2.txt"):
    data = np.loadtxt(file_path, skiprows=1)  # ignorant le titre
    N = data[:,0]
    Erreur_max = data[:,1]

    # Logarithmes des données
    x = np.log(N)
    y = np.log(Erreur_max)

    # Régression linéaire sur les données log-log
    a, b = np.polyfit(x, y, 1)  # y ≈ a*x + b
    beta = -a        # car y = log(ε) = -β*log(N) + log(C)
    C = np.exp(b)

    print(f"Taux de convergence expérimental β ≈ {beta:.3f}")
    print(f"Constante C ≈ {C:.3e}")

    # Tracer les données et la régression linéaire
    plt.figure(figsize=(8,5))
    plt.scatter(x, y, color='blue', label='Données (log(N),log(εN))')
    plt.plot(x, a*x + b, color='red', label=f"Régression linéaire : y = {a:.3f}x + {b:.3f}")
    plt.xlabel("log(N)")
    plt.ylabel("log(εN)")
    plt.title("Convergence expérimentale du Schéma numérique pour l'équation de chaleur 1D")
    plt.grid(True , which="both", ls="--")
    plt.legend()
    plt.show()

plot_erreur("erreur2.txt")