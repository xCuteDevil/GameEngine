# Add sub-directories of your libraries
add_subdirectory("Vector4D")
add_subdirectory("Matrix4x4")
add_subdirectory("Quaternion")
add_subdirectory("PolarCoords")
add_subdirectory("AxisAngle")
add_subdirectory("Vertex")


# Include your liblaries to linking of the project.
target_link_libraries(PA199_project
    PUBLIC 
        lodepng # This links the lodepng library. Do not remove.
        "Vector4D"
        "Matrix4x4"
        "Quaternion"
        "PolarCoords"
        "AxisAngle"
        "Vertex"
)