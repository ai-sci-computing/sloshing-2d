cmake_minimum_required(VERSION 3.20)
project(SloshinTank VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Dependencies ---
find_package(OpenGL REQUIRED)
find_package(glfw3 3.3 REQUIRED)
find_package(doctest REQUIRED)

# --- Main executable ---
add_executable(sloshing_tank
    src/main.cpp
    src/grid.cpp
    src/solver.cpp
    src/vof.cpp
    src/renderer.cpp
    src/simulation.cpp
)

target_include_directories(sloshing_tank PRIVATE src)
target_link_libraries(sloshing_tank PRIVATE OpenGL::GL glfw)

# On macOS, link Cocoa/IOKit frameworks needed by GLFW
if(APPLE)
    target_link_libraries(sloshing_tank PRIVATE
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreVideo"
    )
endif()

# --- Tests ---
enable_testing()

add_executable(test_grid tests/test_grid.cpp src/grid.cpp)
target_include_directories(test_grid PRIVATE src)
target_link_libraries(test_grid PRIVATE doctest::doctest)

add_executable(test_solver tests/test_solver.cpp src/grid.cpp src/solver.cpp src/vof.cpp)
target_include_directories(test_solver PRIVATE src)
target_link_libraries(test_solver PRIVATE doctest::doctest)

add_executable(test_vof tests/test_vof.cpp src/grid.cpp src/vof.cpp)
target_include_directories(test_vof PRIVATE src)
target_link_libraries(test_vof PRIVATE doctest::doctest)

add_executable(test_physics tests/test_physics.cpp src/grid.cpp src/solver.cpp src/vof.cpp)
target_include_directories(test_physics PRIVATE src)
target_link_libraries(test_physics PRIVATE doctest::doctest)

add_test(NAME GridTests COMMAND test_grid)
add_test(NAME SolverTests COMMAND test_solver)
add_test(NAME VofTests COMMAND test_vof)
add_test(NAME PhysicsTests COMMAND test_physics)
