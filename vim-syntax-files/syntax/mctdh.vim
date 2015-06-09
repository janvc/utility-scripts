" Vim syntax file
" Language: MCTDH inputs
" Maintainer: Sebastian Lenz <quizzmaster@gmail.com>
" Latest Revision: 09 Mar 2014

if exists("b:current_syntax")
  finish
endif

" cause using - in keywords
set iskeyword+=-

" Comments
syn keyword mctdhTodo contained TODO FIXME XXX NOTE
syn match mctdhComment "#.*$" contains=mctdhTodo

" Most keywords are case insensitive
syn case ignore

" taken from fortran syntax
" Numbers of various sorts
" Integers
syn match mctdhNumber	display "\<\d\+\(_\a\w*\)\=\>"
" floating point number, without a decimal point
syn match mctdhFloat	display	"\<\d\+[deq][-+]\=\d\+\(_\a\w*\)\=\>"
" floating point number, starting with a decimal point
syn match mctdhFloat	display	"\.\d\+\([deq][-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" floating point number, no digits after decimal
syn match mctdhFloat	display	"\<\d\+\.\([deq][-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" floating point number, E exponents
syn match mctdhFloat	display	"\<\d\+\.\d\+\([E][-+]\=\d\+\)\=\(_\a\w*\)\=\>"
" floating point number
syn match mctdhFloat	display	"\<[-+]\d\+\.\d\+\([deq][-+]\=\d\+\)\=\(_\a\w*\)\=\>"

" Keywords for RUN-SECTION
syn keyword InputRunKeywords contained name stop continue tfinal tout tpsi
syn keyword InputRunKeywords contained psi expect expect1 gridpop steps auto exact
syn keyword InputRunKeywords contained propagation diagonalisation relaxation genoper gendvr geninwf

" Keywords for OPERATOR-SECTION
syn keyword InputOpKeywords contained opname

" Keywords for sbasis-section
syn keyword sbasisKeywords contained multi-set

" Keywords for pbasis-section
syn keyword pbasisKeywords contained FFT HO el rHO Leg Lagu1 Lagu2 Lagu3 Lagu4 sin exp cos sphFBR KLeg K PLeg Wigner Extern xi-xf

" Keywords for INTEGRATOR-SECTION
syn keyword IntegratorKeywords contained VMF CMF ABM BS RK5 RK8 SIL CSIL DAV rDAV rrDAV cDAV

" Keywords for INIT_WF-SECTION
syn keyword InitWFKeywords contained file build end-build gauss HO eigenf pop

" END of Input
syn keyword EndInput end-input end-operator

" operator file sections
" op_define-section
syn keyword OpDefineKeywords contained title end-title

" hamiltonian-section
syn keyword HamiltonianKeywords contained modes KE usediag


" INPUT SECTIONS
syn region runSection matchgroup=mctdhSections start="RUN-SECTION" end="end-run-section" transparent contains=InputRunKeywords,mctdhNumber,mctdhFloat,mctdhComment
syn region OpSection matchgroup=mctdhSections start="OPERATOR-SECTION" end="end-operator-section" transparent contains=InputOpKeywords,mctdhComment
syn region sbasisSection matchgroup=mctdhSections start="sbasis-section" end="end-sbasis-section" transparent contains=mctdhNumber,mctdhFloat,mctdhComment,sbasisKeywords
syn region pbasisSection matchgroup=mctdhSections start="pbasis-section" end="end-pbasis-section" transparent contains=mctdhNumber,mctdhFloat,mctdhComment,pbasisKeywords
syn region IntegratorSection matchgroup=mctdhSections start="INTEGRATOR-SECTION" end="end-integrator-section" transparent contains=mctdhNumber,mctdhFloat,mctdhComment,IntegratorKeywords
syn region InitWFSection matchgroup=mctdhSections start="INIT_WF-SECTION" end="end-init_wf-section" transparent contains=mctdhNumber,mctdhFloat,mctdhComment,InitWFKeywords
" OPERATOR SECTIONS
syn region OpDefineSection matchgroup=mctdhSections start="OP_DEFINE-SECTION" end="end-op_define-section" transparent contains=OpdefineKeywords,mctdhComment
syn region ParameterSection matchgroup=mctdhSections start="PARAMETER-SECTION" end="END-PARAMETER-SECTION" transparent contains=mctdhComment,mctdhNumber,mctdhFloat
syn region HamiltonianSection matchgroup=mctdhSections start="HAMILTONIAN-SECTION" end="END-HAMILTONIAN-SECTION" transparent contains=HamiltonianKeywords,mctdhComment,mctdhNumber,mctdhFloat

" Highlights
hi def link mctdhComment        Comment
hi def link mctdhTodo           Todo
hi def link mctdhNumber         Number
hi def link mctdhFloat          Float
hi def link mctdhFloat          Float
hi def link InputRunKeywords    Type
hi def link InputOpKeywords     Type
hi def link sbasisKeywords      Type
hi def link pbasisKeywords      Type
hi def link IntegratorKeywords  Type
hi def link InitWFKeywords      Type
hi def link OpDefineKeywords    Type
hi def link HamiltonianKeywords Type
hi def link mctdhSections       Statement
hi def link EndInput            Error

syn case match
let b:current_syntax = "mctdh"
