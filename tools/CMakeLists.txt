################################################################################
# Rimerge

link_directories(/usr/local/lib)

add_executable(rimerge.x main.cpp)
add_dependencies(rimerge.x bigbwt)
target_include_directories(rimerge.x PRIVATE ${rimerge-includes})
target_link_libraries(rimerge.x rimerge ${rimerge-libs})

target_compile_features(rimerge.x PUBLIC cxx_constexpr cxx_lambdas cxx_long_long_type cxx_std_17)

################################################################################
# Concatenate FA

add_executable(cfa.x concat_fa.cpp)
target_include_directories(cfa.x PRIVATE ${rimerge-includes})
target_link_libraries(cfa.x rimerge ${rimerge-libs} z)
target_compile_features(cfa.x PUBLIC cxx_constexpr cxx_lambdas cxx_long_long_type cxx_std_17)

################################################################################
# Big-BWT to saes

add_executable(estw.x estw.cpp)
target_include_directories(estw.x PRIVATE ${rimerge-includes})
target_link_libraries(estw.x rimerge ${rimerge-libs})
target_compile_features(estw.x PUBLIC cxx_constexpr cxx_lambdas cxx_long_long_type cxx_std_17)

################################################################################
# Check

add_executable(check.x check_structure.cpp)
target_include_directories(check.x PRIVATE ${rimerge-includes})
target_link_libraries(check.x rimerge ${rimerge-libs})
target_compile_features(check.x PUBLIC cxx_constexpr cxx_lambdas cxx_long_long_type cxx_std_17)

################################################################################
# plain to rle

add_executable(rle.x rle.cpp)
target_include_directories(rle.x PRIVATE ${rimerge-includes})
target_link_libraries(rle.x rimerge ${rimerge-libs})
target_compile_features(rle.x PUBLIC cxx_constexpr cxx_lambdas cxx_long_long_type cxx_std_17)


################################################################################
# stats

add_executable(stats.x stats.cpp)
target_include_directories(stats.x PRIVATE ${rimerge-includes})
target_link_libraries(stats.x rimerge ${rimerge-libs})
target_compile_features(stats.x PUBLIC cxx_constexpr cxx_lambdas cxx_long_long_type cxx_std_17)

################################################################################
# Install
install(TARGETS rimerge.x stats.x rle.x check.x estw.x rimerge.x cfa.x  RUNTIME DESTINATION ${INSTALL_PATH})

