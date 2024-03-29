################################################################################
# Common Framework for Computer Graphics Courses at FI MUNI.
#
# Copyright (c) 2021-2022 Visitlab (https://visitlab.fi.muni.cz)
# All rights reserved.
#
# Course: PA199 (Project Template)
################################################################################

# Generates the lecture. DO NOT MODIFY THIS COMMAND:
visitlab_generate_lecture(PA199 project)

# Finds the external libraries and load their settings.
find_package(lodepng CONFIG REQUIRED)

target_sources(PA199_project
    # Here you can add more source files to the excutable of the game.
    PUBLIC 
    #"some_game_file.cpp" # This includes the tutorial file. Can be removed.
)

# Add sub-directories of your libraries
#add_subdirectory("some_my_library_dir")  # This includes the tutorial library. Can be removed.
#add_subdirectory("Vector4D")
#add_subdirectory("Matrix4x4")
#add_subdirectory("Quaternion")
#add_subdirectory("AxisAngle")

# Include your liblaries to linking of the project.
target_link_libraries(PA199_project
    PUBLIC 
        lodepng # This links the lodepng library. Do not remove.
        #"some_my_library" # This links the tutorial library. Can be removed.
        #"Vector4D"
        #"Matrix4x4"
        #"Quaternion"
        #"AxisAngle"
)

# Includes 'unit_tests' project to the build. DO NOT MODIFY THIS COMMAND:
add_subdirectory("test")  # Do not rename the directory name.
add_subdirectory("Engine")
add_subdirectory("Opengl")

# Installing the application. Do not modify!
install(TARGETS PA199_project DESTINATION .)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/glfw3.dll OPTIONAL DESTINATION .)
install(CODE [[file(WRITE ${CMAKE_INSTALL_PREFIX}/configuration.toml "framework_dir = \".\"\nlecture_dir = \".\"")]])
install(DIRECTORY ./data DESTINATION .)
