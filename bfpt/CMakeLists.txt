message("...................................")
get_filename_component(TARGET_NAME "${CMAKE_CURRENT_LIST_DIR}" NAME)
message("CMake is about to process: ${TARGET_NAME}")

file(GLOB src "src/*.cpp")
message("  Auto detected src files for target ${TARGET_NAME}: ")
foreach(loop_var ${src})
  message("    ${loop_var}")
endforeach(loop_var)
if(NOT src)
  message("    <no files>")
endif()

file(GLOB inc "include/${TARGET_NAME}/*.hpp")
message("  Auto detected hpp files for target ${TARGET_NAME}: ")
foreach(loop_var ${inc})
  message("    ${loop_var}")
endforeach(loop_var)
if(NOT inc)
  message("    <no files>")
endif()

add_executable("${TARGET_NAME}" ${src} ${inc})
target_include_directories("${TARGET_NAME}" PUBLIC include)
target_link_libraries("${TARGET_NAME}"
  Boost::boost
  Boost::program_options
)

message("CMake has processed: ${TARGET_NAME}")