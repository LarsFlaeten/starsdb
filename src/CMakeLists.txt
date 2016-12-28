set(EXENAME starsdb)

# Sources
include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/libraries" "${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB SOURCE_FILES *.cpp)

# Add defines for exe build
# add_definitions("-DXXX=")

link_directories("${PROJECT_SOURCE_DIR}/../cspice/lib")
include_directories("${PROJECT_SOURCE_DIR}/../cspice/include")

message("${PROJECT_SOURCE_DIR}") 

add_executable(${EXENAME} ${SOURCE_FILES})
target_compile_options(${EXENAME} PUBLIC -std=c++11 -Wall)
set(LIBS "cspice")
target_link_libraries(${EXENAME} ${LIBS})

