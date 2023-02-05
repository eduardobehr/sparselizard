set(THIS_APP tracks_electrostatic)

set(CMAKE_BUILD_TYPE Debug)
#find_package(OpenCASCADE REQUIRED)
add_compile_options(-fdiagnostics-color=always)

custom_add_executable_from_dir(${THIS_APP} ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
    ${THIS_APP} sparselizard 
    # gmsh # There are missing dependencies of gmsh
    # occt
    # libocct-foundation-dev
)
custom_symlink_file(${THIS_APP} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} "*.msh;*.geo;*.txt;*.csv;*.nas;*.slz")

