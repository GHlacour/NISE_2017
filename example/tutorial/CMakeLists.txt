# Define executables
add_executable(stochastic EXCLUDE_FROM_ALL src/stochastic.c $<TARGET_OBJECTS:random_lib>)

# Link math
if(NOT WIN32)
    target_link_libraries(stochastic m)
endif()

# Set output directory correct
# We want the compiled code to be put into the current directory for demo purposes
set_target_properties(
    stochastic
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Set aliases
add_executable(Example::Tutorial ALIAS stochastic)