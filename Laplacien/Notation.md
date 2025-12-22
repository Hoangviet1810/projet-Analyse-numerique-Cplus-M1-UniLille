Projet_partie2/
   ├── matrice.h
   ├── matrice.cpp
   ├── matricebande.h
   ├── matricebande.cpp
   ├── main.cpp                 # pour partie 1 & 2
   ├── main2.cpp                # pour partie 3
   └── main_inverse.cpp         # pour partie bonus

├── graph/
    ├── graph.py                # pour partie 1
    ├── graph2.py               # pour partie 3 le cas (1.1)
    ├── graph3.py               # pour partie 3 le cas (1.2)
    └── graph4.py               # pour partie bonus
 
Pour compiler le programme de partie 1 et 2  `main` avec les fichiers `matrice.cpp` et `matricebande.cpp` :

```bash
g++ main.cpp matrice.cpp matricebande.cpp -o test
```

Pour compiler le programme de partie 3  `main2` avec les fichiers `matrice.cpp` et `matricebande.cpp` :

```bash
g++ main2.cpp matrice.cpp matricebande.cpp -o test
```

Pour compiler le programme de partie bonus`main_inverse.cpp` avec les fichiers `matrice.cpp` et `matricebande.cpp` :

```bash
g++ main_inverse.cpp matrice.cpp matricebande.cpp -o test
```

Et ensuite 
```bash
./test
```

Pour compiler le graph
```bash
python3 graph.py
```