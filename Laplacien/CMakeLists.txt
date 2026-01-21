cmake_minimum_required(VERSION 3.5)
project(go)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER g++)

# Fichiers sources communs
set(COMMON_SOURCES
    src/matrice.cpp
    src/matricebande.cpp
)

# Executable 1 : Steepest + Gradient pas fixe
add_executable(go
    src/main.cpp
    ${COMMON_SOURCES}
)

# Executable 2 : PCG
add_executable(go_pcg
    src/main2.cpp
    ${COMMON_SOURCES}
)

# Executable 3 : Probleme inverse
add_executable(go_inverse
    src/main_inverse.cpp
    ${COMMON_SOURCES}
)

# Executable 4 : Test pour tester les fonctions assemblageT, preconditionne et PCG de partie 3
add_executable(go_test
    src/test.cpp
    ${COMMON_SOURCES}
)
# Recherche de l'interpréteur Python3
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Target pour tout exécuter
add_custom_target(run_all
    COMMAND $<TARGET_FILE:go>
    COMMAND $<TARGET_FILE:go_pcg>
    COMMAND $<TARGET_FILE:go_inverse>
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/python/run_graphs.py
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python
    DEPENDS go go_pcg go_inverse go_test
)
