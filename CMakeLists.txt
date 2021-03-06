cmake_minimum_required(VERSION 3.5)

project(grid LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)



find_package(Boost 1.75.0 COMPONENTS date_time timer coroutine thread system unit_test_framework regex REQUIRED)

find_package(Qt5 5.14.2 PATHS /opt/Qt/5.14.2/ COMPONENTS Core Widgets DataVisualization REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

file(GLOB SOURCES
    *.hpp
    *.h
    *.cpp)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -parallel")

add_executable(grid ${SOURCES})

target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX
    Qt5::Core Qt5::Widgets Qt5::DataVisualization ${Boost_LIBRARIES}
)

