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

if ( TARGET_NAME MATCHES "_tests$" )
  message("  Assumed the target is a test target (judged by the target name sufix).")
  string(REGEX REPLACE "_tests" ""
         TESTED_TARGET_NAME "${TARGET_NAME}")
  message("  Tested target: ${TESTED_TARGET_NAME}")
  target_link_libraries("${TARGET_NAME}" ${GTEST_LIBRARIES} pthread)
  target_link_libraries("${TARGET_NAME}" "${TESTED_TARGET_NAME}")
endif()

add_test(NAME "${TARGET_NAME}" COMMAND "${TARGET_NAME}")

message("CMake has processed: ${TARGET_NAME}")