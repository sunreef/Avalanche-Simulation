
cmake_minimum_required(VERSION 2.8)

project(Avalanche)

file(GLOB Sources "src/*.cpp")
file(GLOB Headers "src/*.h")

add_executable(Avalanche ${Sources} ${Headers})

# OpenGL
find_package(OpenGL REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if (UNIX)
  set(GLFW_INCLUDE_DIR /usr/local/Cellar/glfw3)
  set(GLFW_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/GLFW/lib/libglfw3.dylib)
else()
  set(GLFW_INCLUDE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/GLFW/include)
  set(GLFW_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/GLFW/lib/glfw3.lib)
endif()


if (UNIX)
	set(GLEW_INCLUDE_DIR /usr/local/include/)
  set(GLEW_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/GLEW/lib/libGLEW.dylib)
else()
  set(GLEW_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/GLEW/include)
  set(GLEW_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/GLEW/lib/glew32.lib)
endif()

set(GLM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/GLM/include)

include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLFW_INCLUDE_DIR} ${GLEW_INCLUDE_DIR} ${GLM_INCLUDE_DIR})

target_link_libraries(Avalanche ${OPENGL_LIBRARIES} ${GLFW_LIBRARY} ${GLEW_LIBRARY})
