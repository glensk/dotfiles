augroup filetypedetect
au! BufRead,BufNewFile *.rc         setfiletype csh
au! BufRead,BufNewFile makedefs     setfiletype make
au! BufRead,BufNewFile *.biblio     setfiletype biblio
au! BufRead,BufNewFile *.sx         setfiletype sfhingx
au! BufRead,BufNewFile *.F          setfiletype fortran95
au! BufRead,BufNewFile *.f90        setfiletype fortran95
au! BufRead,BufNewFile *.std        setfiletype fhitd
au! BufRead,BufNewFile *.sh         setfiletype sh
au! BufRead,BufNewFile *.math       setfiletype mma
au! BufRead,BufNewFile *.tex        setfiletype tex
augroup END

hi DocComment    guifg=gray70
hi Comment    guifg=gray70

