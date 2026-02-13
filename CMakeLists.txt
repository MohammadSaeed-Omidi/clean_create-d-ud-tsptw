cmake_minimum_required(VERSION 3.10.0)
project(create_tsptw VERSION 0.1.0 LANGUAGES C)

# Find the igraph library
find_library(IGRAPH_LIBRARY igraph REQUIRED)

add_executable(create_tsptw
                main.c
                config.h
                create_tsptw_instances.c
                create_tsptw_instances.h
                solve_instance.h
                solve_instance.c
)

# Link against igraph and math libraries
target_link_libraries(create_tsptw ${IGRAPH_LIBRARY} m)

# Add include directories
target_include_directories(create_tsptw PRIVATE /usr/local/include)
