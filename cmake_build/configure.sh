#!/bin/sh

# this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
# to use it, copy it to another file and run it

# additional compiler flags could be added customizing the corresponding var, for example
# -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 ". Note that the "defaults are already correctly coded"
# so we should add here only machine specific stuff

# an effort is made to autodetect all of the libraries needed HOWEVER:
# METIS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)
# TRILINOS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)

#the user should also note that the symbol "\" marks that the command continues on the next line. IT SHOULD ONLY BE FOLLOWED
#BY the "ENTER" and NOT by any space!!

clear

# you may want to decomment this the first time you compile
rm CMakeCache.txt
rm *.cmake
rm -rf CMakeFiles\

cmake ..                                                                            \
-DCMAKE_C_COMPILER=/usr/bin/gcc                                                     \
-DCMAKE_CXX_COMPILER=/usr/bin/g++                                                   \
-DCMAKE_INSTALL_RPATH="${HOME}/Kratos2019/libs"                                   \
-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE                                            \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 -std=c++11 "                           \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3 "                                          \
-DBOOST_ROOT="${HOME}/boost"                                           \
-DPYTHON_EXECUTABLE="/usr/bin/python3"                                           \
-DCMAKE_BUILD_TYPE=Release                                                       \
-DMESHING_APPLICATION=ON                                                            \
-DEXTERNAL_SOLVERS_APPLICATION=ON                                                   \
-DINCLUDE_FEAST=ON                                                                  \
-DSTRUCTURAL_APPLICATION=OFF                                                         \
-DSTRUCTURAL_MECHANICS_APPLICATION=ON                                               \
-DDEM_STRUCTURES_COUPLING_APPLICATION=OFF                                            \
-DCONVECTION_DIFFUSION_APPLICATION=OFF                                               \
-DSOLID_MECHANICS_APPLICATION=ON                                                    \
-DCONSTITUTIVE_MODELS_APPLICATION=ON                                                \
-DFLUID_DYNAMICS_APPLICATION=OFF                                                     \
-DSHALLOW_WATER_APPLICATION=OFF                                                       \
-DFREE_SURFACE_APPLICATION=OFF                                                       \
-DINCOMPRESSIBLE_FLUID_APPLICATION=OFF                                               \
-DMESH_MOVING_APPLICATION=ON                                                        \
-DMAPPING_APPLICATION=ON                                                           \
-DFSI_APPLICATION=OFF                                                               \
-DMPI_SEARCH_APPLICATION=OFF                                                        \
-DDELAUNAY_MESHING_APPLICATION=ON                                                \
-DCONTACT_MECHANICS_APPLICATION=ON                                               \
-DCONTACT_STRUCTURAL_MECHANICS_APPLICATION=ON                                    \
-DPFEM_APPLICATION=ON                                                           \
-DPFEM_SOLID_MECHANICS_APPLICATION=OFF                                            \
-DPFEM_FLUID_DYNAMICS_APPLICATION=OFF                                              \
-DPFEM2_APPLICATION=OFF                                                           \
-DPARTICLE_MECHANICS_APPLICATION=OFF                                                \
-DLAGRANGIAN_MPM_APPLICATION=OFF                                                    \
-DDEM_APPLICATION=OFF                                                               \
-DSWIMMING_DEM_APPLICATION=OFF                                                      \
-DMIXED_ELEMENT_APPLICATION=OFF                                                     \
-DSHAPE_OPTIMIZATION_APPLICATION=OFF                                                \
-DTOPOLOGY_OPTIMIZATION_APPLICATION=OFF                                             \
-DHDF5_APPLICATION=OFF                                                               \
-DMETIS_APPLICATION=OFF                                                             \
-DPARMETIS_ROOT_DIR="/home/youruser/compiled_libraries/ParMetis-3.1.1"              \
-DTRILINOS_APPLICATION=OFF                                                          \
-DTRILINOS_ROOT="/home/youruser/compiled_libraries/trilinos-10.2.0"                 \
-DINSTALL_EMBEDDED_PYTHON=ON

# decomment this to have it verbose
# make VERBOSE=1 -j4
make -j1
make install
