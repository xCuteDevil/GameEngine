# Add sub-directories of your libraries
add_subdirectory("Camera")
add_subdirectory("Shape")
add_subdirectory("Mesh")

# Include your liblaries to linking of the project.
target_link_libraries(PA199_project
    PUBLIC 
        lodepng # This links the lodepng library. Do not remove.
        "Camera"
        "Shape"
        "Mesh"
)