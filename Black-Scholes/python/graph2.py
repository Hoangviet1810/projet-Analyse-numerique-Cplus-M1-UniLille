# Graph heatmap de l'évolution du payoff du contrat et le prix de l'option d'achat
# pour launcher : il faut de faire commenter la partie 3D et décommenter cette partie

import numpy as np
import matplotlib.pyplot as plt

# Chargement des données
data = np.loadtxt("BS_implicit.txt", comments="#")
T = data[:, 0]        # temps
U = data[:, 1:]       # solution numérique
X = np.linspace(0, 500, U.shape[1])  # espace

# Tracé heatmap X-T
plt.figure(figsize=(10,6))
plt.pcolormesh(X, T, U, shading='auto', cmap='hot')  # Colormap "hot"
cbar = plt.colorbar()
cbar.set_label('Valeur de la solution U(X,T)', fontsize=12)

# Noms des axes
plt.xlabel('Position X (espace)', fontsize=12)
plt.ylabel('Temps T', fontsize=12)

# Titre
plt.title("Évolution du payoff du contrat et le prix de l'option d'achat", fontsize=14)
plt.gca().invert_yaxis()  # Optionnel: T augmente vers le haut si on commente cette ligne
plt.show()