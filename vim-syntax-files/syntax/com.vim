" Vim syntax file
" Language:	Gaussian (Electronic Structure Package) Input Files
" Maintainer:	Jarvist Moore Frost
" Last Change: 21 Aug 2014
" Web Source:	http://jarvistmoorefrost.wordpress.com/2012/05/28/vim-syntax-highlighting-for-nwchem-gaussian-pdb/

" Tweaks:	H.A.Trujillo, 13 Jul 12
" LatestTwkg:	24 Jul 12


if exists("b:current_syntax")
  finish
endif

"Oh Fortran, always ignoring case...
syn case ignore

" Keywords: {{{1

syn keyword	ComBlockCmd dft tddft geometry end basis
syn keyword	ComTask "#.*$"
syn keyword	ComStartup "\%.*$"

syn region	ComDoubleQuote	start=+"+ skip=+\\"+ end=+"+

syn match	ComAtom	"^ \?\a\( \|$\)"

syn match	ComNumber	"-\=\<\d\+\>#\="
syn match	ComNumber	"\<\d\>" display
syn match	ComFloat	"\.\d\+\%([eE][+-]\=\d\+\)\=[jJ]\=\>" display
syn match	ComFloat	"\<\d\+[eE][+-]\=\d\+[jJ]\=\>"  display
syn match	ComFloat	"\<\d\+\.\d*\%([eE][+-]\=\d\+\)\=[jJ]\=" display
"}}}1

" Comments: {{{1
"==========
syn cluster	shCommentGroup  contains=shTodo,@Spell
syn keyword	shTodo          contained	COMBAK FIXME TODO XXX NOTE
syn match	shComment	"^\s*\zs#.*$"   contains=@shCommentGroup
syn match	shComment	"\s\zs#.*$"     contains=@shCommentGroup
syn match	shQuickComment  contained       "#.*$"
syn match 	ComComment	"^\s*!"
"}}}1

" HAT Additions: {{{1
"
"Link-0 fields
sy match	ComChk	"\f\+\.chk"   contains=ComLink0
sy match	ComLink0	"%\a\+="he=e-1	

"Mode Redundant
sy match	ComMR_Type	"^[BADXLO] "	containedin=ComMR_Atoms contained
sy match 	ComMR_Atoms	"^\(\a\s\+\)\?\(\(\d\{1,2}\|\*\)\s\+\)\{2,4}"
"sy match 	ComMR_Atoms	"^\(\a\s\+\)\?\(\(\d\{1,2}\|*\)\s\+\)\{2,4}"

"Geometry Specifications
sy match	ComGeom_Neighb	" \d\{1,2} "	containedin=ComGeom_Rgn	contained
sy match 	ComGeom_Atom	"^ \d\{1,2}"	containedin=ComGeom_Rgn
sy match 	ComGeom_Rgn	"^ \d\{1,2}\(\d\|\s\|\.\)*$"	contains=ComFloat
"}}}1

"Let there be colour {{{1
let b:current_syntax = "com"
hi def link	shTodo	Todo
hi def link	ComBlockCmd	Statement
hi def link	ComTask	Special
hi def link	shDoubleQuote	String
hi def link	shComment	Comment
hi def link	ComNumber	Number
hi def link	ComAtom	PreProc
hi def link	ComStartup	Comment
hi def link	ComFloat	Float

hi link	ComChk	Special
hi link	ComGeom_Atom	PreProc
hi link	ComGeom_Neighb	Statement
"hi link	ComGeom_Neighb	Type
hi link	ComLink0	Statement
hi link	ComMR_Atoms	Statement
hi link	ComMR_Type	Type
hi link ComComment	Comment
"}}}1

" Non_Syntax Settings: {{{1
" ===  Bits that really belong in a ftplugin, not here: === 

" Tweak vim settings
set comments+=nb:!
set textwidth=0

" New Commands
"    :chkpath	to turn filename of chkfile into a full path
"    	       nb/ this assumes that the com-file is being edited from the
"    	           directory containing the chkfile
"
"          1. generic: for use on HOPE only
"cabbrev chkpath     1/chk=/<bar>s/=/=<C-V><C-M>/<bar>-r!pwd<cr>:normal Jr/kJx
"
"          2. tailored for sshfs: 
"                 /Volumes.hope.chem.wilkes.edu/...  --> "                 /home/username/...
"cabbrev chkpath     1/chk=/<bar>s/=/=<C-V><C-M>/<bar>-r!pwd <bar>sed -E 's,/Volumes/hope.(chem<bar>bio).wilkes.edu[^/]*/,/home/username/,'<cr>:normal Jr/kJx<cr>
"
"	   /private/tmp/mydir/...  -->  /home/username/...
cabbrev chkpath     1/chk=/<bar>s/=/=<C-V><C-M>/<bar>-r!pwd<bar>sed s,/private/tmp/mydir,/home/username,<cr>:normal Jr/kJx

"}}}1

" vim: foldmethod=marker ts=15

