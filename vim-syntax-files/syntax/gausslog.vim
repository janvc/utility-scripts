" Vim syntax file
" Language:	Gaussian Log files
" Maintainer:	Jason L. Sonnenberg <sonnenberg.11@osu.edu>
" Last Change:	2011 June 20
" Remark:	Designed to make log files easier to read since 
" 		humans don't generate these files
"
" Acquired Frm:	http://jasonsonnenberg.sufaculty.org/downloads/
" Tweakings:	H.A.Trujillo, 10 Jun 2012
" Last Change:	H.A.Trujillo, 28 Jul 2012

if version < 600
  syntax clear
elseif exists("b:current_syntax")
  finish
endif

" Heading Info	{{{1
"
" At Files 
syntax match glogAtFiles /AtFile.\+/
highlight link glogAtFiles Include

" Copyright
syntax region glogCopyright start=/^ Copyright/ end=/prohibited\sabove/ contains=logDashedLine keepend fold
highlight link glogCopyright Comment

" Citation
syntax region glogCitation start=/^ Cite\sthis\swork/ end=/\.$/	fold
"highlight link glogCitation Type
highlight link glogCitation PreProc	" dark blue

" Link0 commands
syntax match glogLink0 /%\w\+/
highlight link glogLink0 Identifier

" Route line
syntax region glogRoute start=/^\s#/ end=/^\s-\{10,}$/ contains=logDashedLine keepend
highlight link glogRoute Keyword
"}}}1

" Section & Table Demarcation	{{{1
syntax match glogSection /!/
syntax region glogSection start=/-\{4,}/ end=/$/ oneline
syntax region glogSection start=/\*\{4,}/ end=/$/ oneline
syntax region glogSection start=/=\{4,}/ end=/$/ oneline
syntax region glogSection start=/\(Grad\)\{2,}/ end=/$/ oneline
syntax region glogSection start=/\(rotor\)\{2,}/ end=/$/ oneline
syntax region glogSection start=/\(IRC-\)\{2,}/ end=/$/ oneline
highlight link glogSection Comment

" Dashed lines contained in other items
syntax region glogDashedLine start=/-\{2,}/ end=/$/ oneline contained
highlight link glogDashedLine Comment
"}}}1

" Operators	{{{1
syntax match glogOperator /=/
highlight link glogOperator Operator
"}}}1

" Termination	{{{1
"
" Error termination
syntax region glogError start=/^\sError\stermination/ end=/$/
highlight link glogError Error

" Normal termination
syntax region glogSuccess start=/^\sNormal\stermination/ end=/$/
highlight link glogSuccess Type
"}}}1


" ----- HAT Tweaks follow -----


" Heading Info	{{{1
sy match glogHdr 	/^ Input/
sy match glogHdr	/^ Output/
sy match glogFile	"=\f\+"		contains=glogOperator

    hi link glogHdr	Identifier	" green
    hi link glogFile	Special		" purple
"}}}1

" Energies, Convergence {{{1
"
" Convergence
 sy keyword glogConv	YES
 sy keyword glogNoConv	NO
 sy region  glogConvHdr	start="\s*Item\s"	end="Converged?"	

    highlight glogConv		ctermbg=darkgreen	ctermfg=black
    highlight glogNoConv	ctermbg=darkred 	ctermfg=gray
    hi link   glogConvHdr	comment

" Overall Energy
  sy region glogEnergy	start=" Sum of electronic[^=]*Free" end="$" oneline
  sy region glogEnergy	start=" E(\u" end="A.U." oneline
    hi link glogEnergy	statement
"}}}1

" Final Quotation	{{{1
 sy match   glogQtn	/^\s*-\{,2\}\(\u\|[ .,;']\)\+$/	

    hi link glogQtn 	Type		" green
"}}}1

" Transition State Calcs {{{1
"
" Imaginary Frequency 
sy region glogFreqRgn	start="^ Frequencies --" end="$"	oneline	
sy match  glogImagFrq	"-\d\+\.\d\+"	contained containedin=glogSection,glogFreqRgn
sy match  glogImagFrq	"\d\+ imaginary frequencies"	contained
			    \ containedin=glogSection,glogFreqRgn

hi link glogImagFrq special


" IRC
sy match glogIRCsumry	"^\s*Summary of reaction path"
sy match glogIRCsumry	"^\s*Minimum found"

    hi link glogIRCsumry special
"}}}1


" From here to the end should really be in a ft-plugin, not a syntax file:


" Define folds	{{{1
"    			nb/ Folding significantly slows the loading of large files.
set foldmethod=syntax


""			    \ end="\n\s*\(Item\|Error\|\u\{2,} \u\{2}\)"	
syn region glogFold	start="^ GradGradGrad"	
			    \ end="^\s*\(Item\|Error\|^\n\n\)"	
			    \ transparent fold keepend  
			    \ contains=glogConvHdr,glogSection,glogEnergy

syn region glogFold	start="^ IRC-IRC-"
			    \ end="^\(\s*Optimization completed\|\n\n\)"
			    \ transparent fold keepend  
			    \ contains=glogConvHdr,glogSection,glogEnergy

"syn region glogFold	start="^\s\+!\s\+(Angstroms and Degrees)"
syn region glogFold	start="^\s\+!\s\+Optimized Parameters"
			    \ end = "^\s\+\a"me=e-4
			    \ transparent fold keepend  

syn region glogFold	start="^ Test job"	end="^\n\n"
			    \ fold containedin=glodFold

normal zM 
"}}}1

" Mark useful points in file {{{1
"
"    Mark energy line as `e'   ( 'e to return to it )
$
silent! ?^ Sum of electronic\|E(\u?
mark e

"    Mark frequency line as `f'   ( 'f to return to it )
1
silent! /^ Frequencies --/
mark f

"    Place cursor at last convergence check.
$
silent! ?\<\(YES\|NO\)\>?
normal z.
	" mark line as `c'   ( 'c to return here )
mark c
"}}}1

"						       vim: tw=0 foldmethod=marker

