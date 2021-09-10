################################################################################
# Custom Big-BWT
set(BIG-BWT "${CMAKE_SOURCE_DIR}/experiments/Big-Bwt")
add_custom_target(bigbwt COMMAND make WORKING_DIRECTORY ${BIG-BWT})
file(GLOB BIG-BWT_files "${BIG-BWT}/*.x")
set(BIG-BWT_files
        "${BIG-BWT_files}"
        "${BIG-BWT}/bwtparse"
        "${BIG-BWT}/bwtparse64"
        "${BIG-BWT}/remap"
        "${BIG-BWT}/simplebwt"
        "${BIG-BWT}/simplebwt64"
        "${BIG-BWT}/unparse"
        "${BIG-BWT}/unparsz"
        "${BIG-BWT}/pfbwtNT.x"
        "${BIG-BWT}/pfbwtNT64.x"
        ${BIG-BWT_files})
install(PROGRAMS ${BIG-BWT_files} DESTINATION ${INSTALL_PATH})
install(PROGRAMS ${BIG-BWT}/bigbwt DESTINATION ${INSTALL_PATH})

################################################################################
# Custom Online r-index
add_subdirectory(OnlineRlbwt)

################################################################################
# Python testing script
install(PROGRAMS run_tests.py DESTINATION ${INSTALL_PATH})

################################################################################
# PFP++
ExternalProject_Add(pfp
        GIT_REPOSITORY "git@github.com:marco-oliva/pfp.git"
        GIT_TAG "origin/templated_parser"
        CMAKE_ARGS  "-DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR}"
                    "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                    "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
                    "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
                    "-DHTSlib_INCLUDE_DIR=${HTSlib_INCLUDE_DIR}"
                    "-DHTSlib_LIBRARY=${HTSlib_LIBRARY}"
                    "-DCMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH}"
        UPDATE_COMMAND ""
        INSTALL_COMMAND "")
ExternalProject_Get_property(pfp INSTALL_DIR)
install(PROGRAMS  ${INSTALL_DIR}/src/pfp-build/pfp++ DESTINATION ${INSTALL_PATH})
