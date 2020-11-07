message("...................................")
get_filename_component(TARGET_NAME "${CMAKE_CURRENT_LIST_DIR}" NAME)
message("CMake is about to process: ${TARGET_NAME}")

file(GLOB src "src/*.cpp")
message("  Auto detected src files for target ${TARGET_NAME}: ")
message_fancy_list("${src}")

file(GLOB inc "include/${TARGET_NAME}/*.hpp")
message("  Auto detected hpp files for target ${TARGET_NAME}: ")
message_fancy_list("${inc}")

add_library("${TARGET_NAME}" ${src} ${inc})
target_include_directories("${TARGET_NAME}" PUBLIC include)
target_link_libraries("${TARGET_NAME}"
  monostar_hamiltonians
  chainkernel
  so_system
  utility
  armadillo
  Boost::boost
)

set_target_properties("${TARGET_NAME}" PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/")

message("CMake has processed: ${TARGET_NAME}")
