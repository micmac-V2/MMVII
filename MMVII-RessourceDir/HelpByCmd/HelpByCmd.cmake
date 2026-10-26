# Script of the target HelpByCmd, run with  cmake -DSRC_DIR=<this directory> -DWORK_DIR=<scratch directory> [-DPANDOC=<exe>] [-DPDFLATEX=<exe>] -P HelpByCmd.cmake
# (a tool not given is searched in the PATH ; given but empty it is missing).
# One X.pdf per help source, made when it is missing or older than its source, MMVII.sty or Images/ :
#   - X.tex exists                : X.tex -> X.pdf (pdflatex). A hand-written .tex always wins over an X.md.
#   - only X.md exists            : X.md -> X.pdf (pandoc + md2tex.lua + md2tex.template, then pdflatex, with MMVII.sty).
# No .tex is ever written in SRC_DIR from a .md : the intermediate .tex, .aux, .log ... are made in WORK_DIR, the only
# file added to SRC_DIR is the .pdf.
# The files are searched at each run, so adding a .md or a .tex needs no reconfiguration.
# With -DSOFT=ON (target MMVII depends on) a failure is a warning, never an error. Without pdflatex nothing is done,
# without pandoc the help made from a .md is skipped (a warning).
# The PDF are reproducible (fixed date, see below) and not versioned.
cmake_minimum_required(VERSION 3.15)

if(NOT SRC_DIR OR NOT WORK_DIR)
    message(FATAL_ERROR "SRC_DIR and WORK_DIR must be defined")
endif()

# Error, or warning only when SOFT
function(help_fail MSG)
    if(SOFT)
        message(WARNING "${MSG}")
    else()
        message(FATAL_ERROR "${MSG}")
    endif()
endfunction()

if(NOT DEFINED PANDOC)
    find_program(PANDOC pandoc)
endif()
if(NOT DEFINED PDFLATEX)
    find_program(PDFLATEX pdflatex)
endif()
if(NOT PDFLATEX)
    if(NOT SOFT)
        message(FATAL_ERROR "pdflatex not found : help of the commands not generated")
    endif()
    return()   # already warned by the configuration of the main CMakeLists
endif()

set(THIS_SCRIPT "${CMAKE_CURRENT_LIST_FILE}")
set(PANDOC_DEPS "${SRC_DIR}/md2tex.lua" "${SRC_DIR}/md2tex.template" "${THIS_SCRIPT}")
file(MAKE_DIRECTORY "${WORK_DIR}")

# Names of the help : every .tex, and every .md without .tex
file(GLOB TEX_FILES "${SRC_DIR}/*.tex")
file(GLOB MD_FILES "${SRC_DIR}/*.md")
file(GLOB IMAGE_FILES "${SRC_DIR}/Images/*")
set(NAMES "")
foreach(F ${TEX_FILES})
    get_filename_component(NAME "${F}" NAME_WLE)
    list(APPEND NAMES "${NAME}")
endforeach()
foreach(F ${MD_FILES})
    get_filename_component(NAME "${F}" NAME_WLE)
    if(NOT EXISTS "${SRC_DIR}/${NAME}.tex")
        list(APPEND NAMES "${NAME}")
    endif()
endforeach()
list(SORT NAMES)

if(WIN32)
    set(PATH_SEP ";")
else()
    set(PATH_SEP ":")
endif()

# Same source, same PDF : fixed date (2025-01-01) instead of the time of the run, and no /ID (pdfTeX derives it
# from the output path, which differs from one build directory to another)
set(PDF_DATE 1735689600)

set(FAILED "")
foreach(NAME ${NAMES})
    set(PDF "${SRC_DIR}/${NAME}.pdf")
    set(TMP "${WORK_DIR}/${NAME}")

    if(EXISTS "${SRC_DIR}/${NAME}.tex")
        set(SRC "${SRC_DIR}/${NAME}.tex")
        set(FROM_MD FALSE)
        set(DEPS "${SRC}")
    else()
        set(SRC "${SRC_DIR}/${NAME}.md")
        set(FROM_MD TRUE)
        set(DEPS "${SRC}" ${PANDOC_DEPS})
    endif()

    if(EXISTS "${PDF}")
        set(UP_TO_DATE TRUE)
        foreach(DEP ${DEPS} "${SRC_DIR}/MMVII.sty" ${IMAGE_FILES})
            if("${DEP}" IS_NEWER_THAN "${PDF}")
                set(UP_TO_DATE FALSE)
            endif()
        endforeach()
        if(UP_TO_DATE)
            continue()
        endif()
    endif()

    file(REMOVE_RECURSE "${TMP}")
    file(MAKE_DIRECTORY "${TMP}")

    # The file given to pdflatex : the .tex itself, or the one pandoc makes in TMP
    if(FROM_MD)
        if(NOT PANDOC)
            help_fail("pandoc not found : ${NAME}.pdf not made from ${NAME}.md")
            continue()
        endif()
        message(STATUS "${NAME}.md -> ${NAME}.pdf")
        execute_process(
            COMMAND "${PANDOC}" -f markdown-auto_identifiers -t latex --wrap=preserve -s
                    "--lua-filter=${SRC_DIR}/md2tex.lua" "--template=${SRC_DIR}/md2tex.template" -V cmd=${NAME}
                    -o "${TMP}/${NAME}.tex" "${SRC}"
            RESULT_VARIABLE RES ERROR_VARIABLE ERR)
        if(NOT RES EQUAL 0)
            help_fail("pandoc failed on ${NAME}.md :\n${ERR}")
            list(APPEND FAILED "${NAME}")
            file(REMOVE_RECURSE "${TMP}")
            continue()
        endif()
        set(TEX_INPUT "${TMP}/${NAME}.tex")
    else()
        message(STATUS "${NAME}.tex -> ${NAME}.pdf")
        set(TEX_INPUT "${NAME}.tex")
    endif()

    set(RES 0)
    foreach(PASS 1 2)   # second pass : table of contents and references
        if(RES EQUAL 0)
            # run in SRC_DIR (relative paths of images and MMVII.sty), output elsewhere ; the trailing separator keeps the default TeX paths
            execute_process(
                COMMAND "${CMAKE_COMMAND}" -E env "TEXINPUTS=${SRC_DIR}${PATH_SEP}" "SOURCE_DATE_EPOCH=${PDF_DATE}" FORCE_SOURCE_DATE=1
                        "${PDFLATEX}" -interaction=batchmode -halt-on-error "-jobname=${NAME}" "-output-directory=${TMP}"
                        "\\pdftrailerid{}\\input{${TEX_INPUT}}"
                WORKING_DIRECTORY "${SRC_DIR}"
                RESULT_VARIABLE RES OUTPUT_QUIET ERROR_QUIET)
        endif()
    endforeach()
    if(RES EQUAL 0)
        execute_process(COMMAND "${CMAKE_COMMAND}" -E copy "${TMP}/${NAME}.pdf" "${PDF}")
    else()
        message(WARNING "pdflatex failed on ${NAME}, see the first error below")
        if(EXISTS "${TMP}/${NAME}.log")
            execute_process(COMMAND grep -m 1 -A 4 "^!" "${TMP}/${NAME}.log")
        endif()
        list(APPEND FAILED "${NAME}")
    endif()
    file(REMOVE_RECURSE "${TMP}")
endforeach()

file(REMOVE_RECURSE "${WORK_DIR}")
if(FAILED)
    help_fail("Not compiled : ${FAILED}")
endif()
