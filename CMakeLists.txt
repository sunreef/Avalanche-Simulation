
cmake_minimum_required(VERSION 2.8)

project(Avalanche)

file(GLOB Sources "src/*.cpp")
file(GLOB Headers "src/*.h")

add_executable(Avalanche ${Sources} ${Headers})

# OpenGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

set(GLFW_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/GLFW/include)
if (UNIX)
  set(FREEGLUT_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/FreeGLUT/lib/libglut.dylib) 
else()
  set(GLFW_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/GLFW/lib/glfw3.lib)
endif()

include_directories(${GLUT_INCLUDE_DIRS})

set(GLEW_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/GLEW/include)
if (UNIX)
  set(GLEW_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/GLEW/lib/libGLEW.dylib)
else()
  set(GLEW_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/GLEW/lib/glew32.lib)
endif()

set(GLM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/GLM/include)

include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLFW_INCLUDE_DIR} ${GLEW_INCLUDE_DIR} ${GLM_INCLUDE_DIR})

target_link_libraries(Avalanche ${OPENGL_LIBRARIES} ${GLFW_LIBRARY} ${GLEW_LIBRARY})
