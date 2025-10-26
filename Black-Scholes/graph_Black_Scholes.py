# Graph en 3D de l'évolution du payoff du contrat et le prix de l'option d'achat
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Chargement des données
data = np.loadtxt("BS_implicit.txt", comments="#")
T = data[:, 0]        # temps
U = data[:, 1:]       # solution numérique
X = np.linspace(0, 500, U.shape[1])  # espace

# Création de la grille X-T
X_grid, T_grid = np.meshgrid(X, T)

# Tracé de la surface 3D
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X_grid, T_grid, U, cmap='hot', edgecolor='none')

# Colorbar indiquant la valeur de la solution
cbar = fig.colorbar(surf, ax=ax, shrink=0.6, aspect=10)
cbar.set_label('Valeur de la solution U(X,T)', fontsize=12)

# Noms des axes
ax.set_xlabel('Position X (espace)', fontsize=12)
ax.set_ylabel('Temps T', fontsize=12)
ax.set_zlabel('Hauteur = Valeur U(X,T)', fontsize=12)

# Titre
ax.set_title("Évolution du payoff du contrat et le prix de l'option d'achat", fontsize=14)

# Vue 3D
ax.view_init(elev=30, azim=-60)

plt.show()



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