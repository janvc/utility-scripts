" Vim syntax file
" Language: Gaussian09 input files (*.com)
" Maintainer: Jan von Cosel <janvc@gmx.de>
" Latest Revision: 80 Jun 2015

if exists("b:current_syntax")
  finish
endif

"Oh Fortran, always ignoring case...
syn case ignore

" taken from fortran syntax
" Numbers of various sorts
" Integers
syn match gaussianNumber	display "\<\d\+\(_\a\w*\)\=\>"
" floating point number, without a decimal point
syn match gaussianFloat	    display	"\<\d\+[deq][-+]\=\d\+\(_\a\w*\)\=\>"
" floating point number, starting with a decimal point
syn match gaussianFloat	    display	"\.\d\+\([deq][-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" floating point number, no digits after decimal
syn match gaussianFloat	    display	"\<\d\+\.\([deq][-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" floating point number, E exponents
syn match gaussianFloat	    display	"\<\d\+\.\d\+\([E][-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" floating point number
syn match gaussianFloat	    display	"\<[-+]\d\+\.\d\+\([deq][-+]\=\d\+\)\=\(_\a\w*\)\=\>"


" Keywords for Link0 section
syn keyword Link0Keywords chk= nprocshared= mem=
syn match Link0Keywords "%\a\+="

" Keywords for Route section:
" output verbosity setting
syn keyword RouteKeywords p t n
" methods
syn keyword RouteKeywords hf mp2 mp3 mp4 wb97xd b3lyp
" basis sets
syn keyword RouteKeywords gen sto-3g svp tzvp
" keyword options: opt
syn keyword RouteKeywords opt tight maxstep
" keyword options: freq
syn keyword RouteKeywords freq hpmodes noraman savenormalmodes savenm
" keyword options: other keywords
syn keyword RouteKeywords formcheck


" sections
syn region RouteSection matchgroup=RouteSections start="#" end="^$" contains=RouteKeywords
syn region Link0Section matchgroup=Link0Sections start="%" end="$" contains=Link0Keywords

hi def link RouteSections  Type
hi def link Link0Sections  Statement
hi def link RouteKeywords  Type
hi def link Link0Keywords  Statement
hi def link gaussianNumber Number
hi def link gaussianFloat  Float
