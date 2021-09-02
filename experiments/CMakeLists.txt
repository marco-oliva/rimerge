################################################################################
# Custom Big-BWT
set(BIG-BWT "${CMAKE_SOURCE_DIR}/experiments/Big-Bwt")
add_custom_target(bigbwt COMMAND make  WORKING_DIRECTORY ${BIG-BWT})
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
        "${BIG-BWT}/pscan.x"
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