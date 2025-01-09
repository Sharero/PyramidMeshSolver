cmake_minimum_required(VERSION 3.10)

project(test)

add_executable(my_app src/main.cpp src/grid.cpp src/fem.cpp src/slae.cpp src/msg.cpp)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_custom_target(visualize
    COMMAND ${CMAKE_COMMAND} -E env python3 ${CMAKE_SOURCE_DIR}/python/draw.py
    DEPENDS my_app
)

# Связываем запуск программы и визуализацию
add_custom_command(
    TARGET my_app POST_BUILD
    COMMAND my_app
    # COMMAND ${CMAKE_COMMAND} --build . --target visualize
    COMMAND ${CMAKE_COMMAND} -E remove $<TARGET_FILE:my_app>
)
