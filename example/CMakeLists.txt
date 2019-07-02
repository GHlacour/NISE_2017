project(NISE_SAMPLES VERSION 1.0 LANGUAGES C)
add_subdirectory("tutorial")
add_subdirectory("tutorialSFG")

add_custom_target(examples)
add_dependencies(examples Example::Tutorial Example::TutorialSFG)