cmake_minimum_required(VERSION 3.14)

project(swe_simulation)

include(FetchContent)

FetchContent_Declare(
    glew
    GIT_REPOSITORY https://github.com/Perlmint/glew-cmake.git
    GIT_TAG        f456deace7b408655109aaeff71421ef2d3858c6 # 2.2.0
)
set(glew-cmake_BUILD_SHARED OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glew)
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        b0796109629931b6fa6e449c15a177845256a407 # 3.3
)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)


# imgui
add_library(imgui STATIC
  imgui/imgui.cpp
  # imgui/imgui_demo.cpp
  imgui/imgui_draw.cpp
  imgui/imgui_widgets.cpp
  imgui/imgui_tables.cpp
  imgui/imgui_impl_glfw.cpp
  imgui/imgui_impl_opengl3.cpp
)

target_include_directories(imgui
PRIVATE 
  ${glew_SOURCE_DIR}/include
)
target_link_libraries(imgui
  libglew_static
  glfw
)

add_executable(${PROJECT_NAME} 
  shallow_water.cpp
)

target_link_libraries(${PROJECT_NAME} 
  imgui
)
