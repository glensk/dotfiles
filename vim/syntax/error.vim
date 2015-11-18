syntax region  errLine        start="lin" end=":"

syntax region  errString       start=+"+ end=+"+

syntax region  errErrInt       start="err" end=":"

syntax region  errFile         start=+"+ end=+"+

if !exists("did_error_syntax_inits")
  let did_error_syntax_inits = 1
  highlight errErrInt	term=underline ctermfg=6 guifg=grey50 gui=bold
  highlight errString  	term=underline ctermfg=6 guifg=blue4 gui=bold
  highlight errLine   	term=underline ctermfg=6 guifg=red4 gui=bold
  highlight errLineNo	term=underline ctermfg=6 guifg=brown gui=bold
endif
