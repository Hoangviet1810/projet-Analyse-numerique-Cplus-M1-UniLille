# Nom de l'exécutable
EXEC = go

# Liste des sources
SRCS = src/main.cpp src/fun.cpp

# Liste des en-têtes
HDRS = src/fun.h

# Compilation et édition des liens
$(EXEC): main.o fun.o
	$(CXX) -o $(EXEC) main.o fun.o

# Règle pour main.o (dépend de main.cpp et fun.h)
main.o: src/main.cpp src/fun.h
	$(CXX) -c src/main.cpp

# Règle pour fun.o (dépend de fun.cpp et fun.h)
fun.o: src/fun.cpp src/fun.h
	$(CXX) -c src/fun.cpp

# Option de nettoyage
clean:
	rm -f *.o $(EXEC)