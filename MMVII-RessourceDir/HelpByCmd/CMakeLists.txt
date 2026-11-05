# Targets making the .pdf of each help of this directory, from its .tex or else from its .md (see HelpByCmd.cmake).
# Included by the main CMakeLists on Unix only, which also finds pandoc and pdflatex (HELP_PANDOC, HELP_PDFLATEX).
#   HelpByCmd      : ninja HelpByCmd ; fails when a tool is missing or a file does not compile
#   HelpByCmd_opt  : same, but only warns ; MMVII depends on it, so any build of MMVII updates the help
foreach(SUFFIX "" "_opt")
    if(SUFFIX STREQUAL "")
        set(SOFT_OPT OFF)
    else()
        set(SOFT_OPT ON)
    endif()
    add_custom_target(HelpByCmd${SUFFIX}
        COMMAND "${CMAKE_COMMAND}" "-DSRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}" "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/tmp${SUFFIX}"
                "-DSOFT=${SOFT_OPT}" "-DPANDOC=${HELP_PANDOC}" "-DPDFLATEX=${HELP_PDFLATEX}" -P "${CMAKE_CURRENT_SOURCE_DIR}/HelpByCmd.cmake"
        COMMENT "Help of the commands : .tex or .md -> .pdf"
        VERBATIM USES_TERMINAL)
endforeach()
