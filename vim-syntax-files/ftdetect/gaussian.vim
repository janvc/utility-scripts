fun! s:DetectGausslog()
    if getline(1) =~ 'Entering Gaussian System.*'
        set filetype=gausslog
    endif
endfun

au BufRead,BufNewFile *.com set filetype=gausscom
au BufRead,BufNewFile *.log call s:DetectGausslog()
