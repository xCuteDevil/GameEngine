# This is the executable into which you are
# supposed to include tests of your math
# library functions.

find_package(GTest CONFIG REQUIRED)

add_executable(PA199_project_tests
    "some_test.cpp" # This includes the tutorial file. Can be removed.
    )

target_link_libraries(PA199_project_tests PUBLIC
  GTest::gtest
  GTest::gtest_main
  # Here add libraries you want to write tests for.
  #some_my_library # This includes the tutorial library. Can be removed.
  Vector4D
  Matrix4x4
  Quaternion
  AxisAngle
  PolarCoords
)

# Exposes the test (e.g., for Visual Studio).
add_test(NAME PA199_project_tests COMMAND PA199_project_tests)

# Installing the tests. Do not modify!
install(TARGETS PA199_project_tests DESTINATION .)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gtest.dll OPTIONAL DESTINATION .)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gtest_main.dll OPTIONAL DESTINATION .)
