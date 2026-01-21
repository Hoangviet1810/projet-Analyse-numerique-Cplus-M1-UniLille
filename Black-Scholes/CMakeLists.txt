cmake_minimum_required(VERSION 3.5)
project(go)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER g++)

# Fichiers sources communs
set(COMMON_SOURCES
    src/Black_Scholes.cpp
    src/Black_Scholes.h
)

# Executable 1 : Steepest + Gradient pas fixe
add_executable(go
    src/main.cpp
    ${COMMON_SOURCES}
)

# Recherche de l'interpréteur Python3
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Target pour tout exécuter
add_custom_target(run_all
    COMMAND $<TARGET_FILE:go>
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/python/run_graphs.py
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python
    DEPENDS go 
)
