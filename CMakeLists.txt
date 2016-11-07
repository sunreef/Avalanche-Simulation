
cmake_minimum_required(VERSION 2.8)

project(Avalanche)

file(GLOB Sources "src/*.cpp")
file(GLOB Headers "src/*.h")

add_executable(Avalanche ${Sources} ${Headers})
find_package(OpenGL REQUIRED)

set(FREEGLUT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/FreeGLUT/include)
set(FREEGLUT_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/FreeGLUT/lib/freeglut.lib ${CMAKE_CURRENT_SOURCE_DIR}/FreeGLUT/lib/freeglutd.lib)

set(GLEW_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/GLEW/include)
set(GLEW_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/GLEW/lib/glew32.lib)

set(GLM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/GLM/include)

include_directories( ${OPENGL_INCLUDE_DIRS}  ${FREEGLUT_INCLUDE_DIR} ${GLEW_INCLUDE_DIR} ${GLM_INCLUDE_DIR})

target_link_libraries(Avalanche ${OPENGL_LIBRARIES} ${FREEGLUT_LIBRARY} ${GLEW_LIBRARY})