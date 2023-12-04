# This is a tutorial file. Feel free to remove it.

add_library(Shape
    
    "Shape.hpp"
    "Shape.cpp"
)

find_package(glad CONFIG REQUIRED)
target_link_libraries(Shape PUBLIC
glad::glad)
