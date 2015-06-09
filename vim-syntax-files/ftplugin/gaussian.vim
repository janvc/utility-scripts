" Vim settings file
" Language:	Gaussian files
" Version:	0.1
" Last Change:	2015 Jun. 09
" Maintainer:	Sebastian Lenz <sebastian@slenz.net>
" Usage:	Do :help gaussian-plugin from Vim
" Credits:

" Only do these settings when not done yet for this buffer
if exists("b:did_ftplugin")
  finish
endif

" Don't do other file type settings for this buffer
let b:did_ftplugin = 1

let b:gaussian_firstline = getline(1)

if b:gaussian_firstline =~# "^Entering Gaussian System, Link 0"
    au BufRead,BufNewFile *.log set filetype=gausslog
endif
