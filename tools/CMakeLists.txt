################################################################################
# Temporary main

link_directories(/usr/local/lib)

add_executable(main main.cpp)
add_dependencies(main rimerge CLI11 mio sdsl)
target_include_directories(main PRIVATE ${rimerge-includes})
target_link_libraries(main rimerge ${rimerge-libs})

target_compile_features(main PUBLIC cxx_constexpr cxx_lambdas cxx_long_long_type cxx_std_17)

################################################################################
# Concatenate FA

add_executable(cfa concat_fa.cpp)
target_include_directories(cfa PRIVATE ${rimerge-includes})
target_link_libraries(cfa rimerge ${rimerge-libs})
target_compile_features(cfa PUBLIC cxx_constexpr cxx_lambdas cxx_long_long_type cxx_std_17)

################################################################################
# HiperGator specific install

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    install(TARGETS main cfa CONFIGURATIONS Release RUNTIME DESTINATION /ufrc/boucher/marco.oliva/data/r-index/tools)
    add_custom_command(TARGET main POST_BUILD COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install)
    add_custom_command(TARGET cfa POST_BUILD COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install)
endif()
