find_package(LATEX)
if(NOT DEFINED GENERATE_DOCS)
    if($<CONFIG:Debug>)
        set(GENERATE_DOCS $<BOOL:false>)
    else()
        set(GENERATE_DOCS $<BOOL:true>)
    endif()
endif()

if(GENERATE_DOCS AND LATEX_COMPILER AND BIBTEX_COMPILER)
    include(UseLATEX)

    add_latex_document(
        NISE.tex
        INPUTS tableofcontents.tex developers.tex techniques.tex input.tex
        IMAGES 2ds.pdf
        BIBFILES Merged2016.bib
        TARGET_NAME doc
        EXCLUDE_FROM_ALL
        EXCLUDE_FROM_DEFAULTS
    )
else()
    message(WARNING "Docs cannot be compiled since no working LaTeX distribution could be found!")
endif()
