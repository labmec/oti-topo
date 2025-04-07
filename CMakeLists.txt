cmake_minimum_required(VERSION 3.14)

project (OtiTopo LANGUAGES C CXX)


set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)

# Making pz "seeable" in Xcode
file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
     "../NeoPzDevelop/*.h"
     "../NeoPzDevelop/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})

# include_directories(
# 	"${PROJECT_SOURCE_DIR}"
# )


# creating the projects
add_executable(oti-topo otitopo.cpp TPZElasticity2DOtiTopoMP.cpp TPZElasticity2DOtiTopo.cpp TPZElasticity3DOtiTopo.cpp TPZHelmholtz.cpp TPZHelmholtz.cpp)
target_link_libraries(oti-topo PUBLIC NeoPZ::pz)

add_executable(oti-topo-matwithmem otitopo_matwithmem.cpp TPZElasticity2DOtiTopo.cpp TPZElasticity3DOtiTopo.cpp TPZHelmholtz.cpp)
target_link_libraries(oti-topo-matwithmem PUBLIC NeoPZ::pz)

add_executable(darcyzinho darcyzinho.cpp TPZElasticity2DOtiTopo.cpp TPZElasticity3DOtiTopo.cpp TPZHelmholtz.cpp)
target_link_libraries(darcyzinho PUBLIC NeoPZ::pz)


