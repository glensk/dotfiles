syn clear

syn case ignore

"if exists("fortran_no_tab_highlight")
  syn match fortTab	"\t" transparent
"else
"  syn match fortTab	"\t"
"endif

syn match fortWS	contained "\s\+" contains=fortTab,fortRightMargin

syn keyword fortTodo 	contained TODO

syn match fortNumber	"\<[0-9]\+\>"	contains=fortRightMargin
syn match fortFloat	"\<[0-9]\+\.[0-9]*\([EDQ][-+]\=[0-9]\+\)\=\>" contains=fortRightMargin
syn match fortFloat	"\.[0-9]\+\([EDQ][-+]\=[0-9]\+\)\=\>" contains=fortRightMargin
syn match fortFloat	"\<[0-9]\+[EDQ][-+]\=[0-9]\+\>" contains=fortRightMargin

syn match fortOperator	"\.\(gt\|ge\|\lt\|le\|eq\|ne\)\." contains=fortRightMargin
syn match fortOperator	"\.\(eqv\|neqv\|and\|or\|not\)\." contains=fortRightMargin
syn match fortOperator	"[-+=><:%]"

syn match fortLogical	"\.\(true\|false\)\." contains=fortRightMargin

"syn match fortIdentifier	"\<[a-z_][a-z0-9_]*\>"	contains=fortRightMargin

syn match fortDelimiter	"[,;]"

" END cannot be a keyword as it is used in matches later on
syn match   fortUnitHeader	"\<end\>" contains=fortRightMargin
syn keyword fortUnitHeader	call entry result
syn match   fortUnitHeader	"\<\(end\s*\)\=block\s*data\>" contains=fortTab,fortRightMargin
syn match   fortUnitHeader	"\<\(end\s*\)\=\(function\|module\|program\|subroutine\)\>" contains=fortTab,fortRightMargin
syn match   fortUnitHeader	"\<contains\>" contains=fortRightMargin

syn keyword fortStatement	continue pause return stop

syn keyword fortRepeat		cycle exit while
syn match   fortRepeat		"\<\(end\s*\)\=\(do\|forall\)\>" contains=fortTab,fortRightMargin

syn match   fortGoTo		"\<go\s*to\>" contains=fortTab,fortRightMargin
syn keyword fortGoTo		assign to

syn keyword fortType		byte character complex integer logical real
syn match   fortType		"\<double\s*\(complex\|precision\)\>" contains=fortTab,fortRightMargin

" Fortran 77 declarations
syn keyword fortDeclaration	common data dimension equivalence external
syn keyword fortDeclaration	implicit intrinsic parameter save

" Fortran 90/95 declarations and attributes
syn keyword fortDeclaration	allocatable elemental in inout intent namelist
syn keyword fortDeclaration	none only operator optional out pointer private
syn keyword fortDeclaration	public pure recursive sequence target use
syn match   fortDeclaration	"\<\(end\s*\)\=\(interface\|type\)\>" contains=fortTab,fortRightMargin
syn match   fortDeclaration	"\<module\s\+procedure\>"	contains=fortTab,fortRightMargin

syn keyword fortConditional	elsewhere if then
syn match   fortConditional	"\<else\>"			contains=fortRightMargin
syn match   fortConditional	"\<\(else\|end\)\s*if\>"	contains=fortTab,fortRightMargin
syn match   fortConditional	"\<\(end\s*\)\=select\>"	contains=fortTab,fortRightMargin
syn match   fortConditional	"\<case\(\s\+default\)\=\>"	contains=fortTab,fortRightMargin
syn match   fortConditional	"\<select\s\+case\>"		contains=fortTab,fortRightMargin
syn match   fortConditional	"\<\(end\s*\)\=where\>"		contains=fortTab,fortRightMargin

syn keyword fortIOControlSpec	contained access action advance blank delim direct end
syn keyword fortIOControlSpec	contained eor err exist file fmt form formatted
syn keyword fortIOControlSpec	contained iolength iostat name named nextrec nml
syn keyword fortIOControlSpec	contained number opened pad position read readwrite rec
syn keyword fortIOControlSpec	contained recl sequential size status unformatted unit write

syn region fortFormat		contained matchgroup=fortDelimiter start="(" end=")" contains=fortTab,fortRightMargin,fortFormat,fortHollerith,fortString,fortDelimiter,fortContMark,fortContError,fortContComment

" the following groups are special to fortIO
syn match fortIOListSkip	contained "!.*$"		skipempty nextgroup=fortIOListSkip contains=fortComment
syn match fortIOListSkip	contained "^\([c*]\|\s*!\).*$"	skipempty nextgroup=fortIOListSkip contains=fortComment,fortTab
syn match fortIOListSkip	contained "^ \{5}.\s*"	nextgroup=fortIOList contains=fortContError,fortContMark

syn match fortFormatSkip	contained "!.*$"		skipempty nextgroup=fortFormatSkip contains=fortComment
syn match fortFormatSkip	contained "^\([c*]\|\s*!\).*$"	skipempty nextgroup=fortFormatSkip contains=fortComment,fortTab
syn match fortFormatSkip	contained "^ \{5}.\s*"	nextgroup=fortFormat contains=fortContError,fortContMark

syn match fortFmtStringSkip	contained "!.*$"		skipempty nextgroup=fortFmtStringSkip contains=fortComment
syn match fortFmtStringSkip	contained "^\([c*]\|\s*!\).*$"	skipempty nextgroup=fortFmtStringSkip contains=fortComment,fortTab
syn match fortFmtStringSkip	contained "^ \{5}.\s*"	nextgroup=fortFmtString contains=fortContError,fortContMark

syn match fortReadSkip		contained "!.*$"		skipempty nextgroup=fortReadSkip contains=fortComment
syn match fortReadSkip		contained "^\([c*]\|\s*!\).*$"	skipempty nextgroup=fortReadSkip contains=fortComment,fortTab
syn match fortReadSkip		contained "^ \{5}.\s*"	nextgroup=fortRead contains=fortContError,fortContMark

syn match fortIO	"\<\(backspace\|close\|inquire\|open\|rewind\|write\)\>\s*" skipempty nextgroup=fortIOList,fortIOListSkip contains=fortWS,fortRightMargin
syn match fortIO	"\<format\>\s*"		skipempty nextgroup=fortFormat,fortFormatSkip contains=fortWS,fortRightMargin
syn match fortIO	"\<end\s*file\>\s*"	skipempty nextgroup=fortIOList,fortIOListSkip contains=fortWS,fortRightMargin
syn match fortIO	"\<print\>\s*"		skipempty nextgroup=fortFmtString,fortFmtStringSkip contains=fortWS,fortRightMargin
syn match fortIO	"\<read\>\s*"		skipempty nextgroup=fortIOList,fortFmtString,fortReadSkip contains=fortWS,fortRightMargin

syn keyword fortIntrinsic	abs acos aint atan asin cos cosh aimag anint
syn keyword fortIntrinsic	atan2 char cmplx conjg dble dim dprod exp
syn keyword fortIntrinsic	ichar index int len lge lgt lle llt log log10
syn keyword fortIntrinsic	max min mod nint sin sinh sign sqrt tan tanh

" Fortran 90/95 intrinsic functions
syn keyword fortIntrinsic	achar adjustl adjustr all allocated any
syn keyword fortIntrinsic	associated bit_size btest ceiling count cpu_time
syn keyword fortIntrinsic	cshift date_and_time digits dot_product eoshift
syn keyword fortIntrinsic	epsilon exponent floor fraction huge iachar
syn keyword fortIntrinsic	iand ibclr ibits ibset ieor ior ishft ishftc
syn keyword fortIntrinsic	kind lbound len_trim matmul maxexponent maxloc
syn keyword fortIntrinsic	maxval merge minexponent minloc minval modulo
syn keyword fortIntrinsic	mvbits nearest not null pack precision present
syn keyword fortIntrinsic	product radix random_number random_seed range
syn keyword fortIntrinsic	repeat reshape rrspacing scale scan
syn keyword fortIntrinsic	selected_int_kind selected_real_kind
syn keyword fortIntrinsic	set_exponent shape size spacing spread sum
syn keyword fortIntrinsic	system_clock tiny transfer transpose trim
syn keyword fortIntrinsic	ubound unpack verify

" Fortran 77 specific functions
syn keyword fortSpecific	alog alog10 amax0 amax1 amin0 amin1 amod cabs
syn keyword fortSpecific	ccos cexp clog csin csqrt dabs dacos dasin
syn keyword fortSpecific	datan datan2 dcos dcosh ddim dexp dint dlog
syn keyword fortSpecific	dlog10 dmax1 dmin1 dmod dnint dsign dsin
syn keyword fortSpecific	dsinh dsqrt dtan dtanh float iabs idim idint
syn keyword fortSpecific	idnint ifix isign max0 max1 min0 min1 sngl

" Some other Fortran 90 syntax
syn keyword fortMemory		allocate deallocate nullify stat
syn keyword fortInclude		include

" Not in any Fortran standard
syn keyword fortExtended	carg carriagecontrol flush
syn keyword fortExtended	getcl system val recordtype

syn keyword fortExtended	blocksize break carg dvchk error from intrup
syn keyword fortExtended	invalop iostat_msg location nbreak ndperr ndpexc
syn keyword fortExtended	offset ovefl precfill prompt segment timer undfl

syn match fortComment	        "!.*$"	
syn match fortDocComment	"!!.*$"	

syn region fortString	start=+"+ end=+"+ contains=fortTab,fortRightMargin,fortContMarkString,fortContError,fortContComment
syn region fortString	start=+'+ end=+'+ contains=fortTab,fortRightMargin,fortContMarkString,fortContError,fortContComment

" Format strings are used with READ, WRITE and PRINT
syn region  fortFmtString	contained start=+'(+ end=+)'+ contains=fortTab,fortRightMargin,fortContMarkString,fortContError,fortContComment
syn region  fortFmtString	contained start=+"(+ end=+)"+ contains=fortTab,fortRightMargin,fortContMarkString,fortContError,fortContComment

syn match  fortParenError	")"
syn region fortParen		matchgroup=fortDelimiter start="(" end=")" transparent contains=fortAmpersandError,fortComment,fortContComment,fortContError,fortContMark,fortDelimiter,fortExtended,fortFloat,fortHollerith,fortIntrinsic,fortNumber,fortOperator,fortParen,fortRightmargin,fortSpecific,fortString,fortTab,fortLogical
syn region fortIOList		contained matchgroup=fortDelimiter start="(" end=")" transparent contains=fortAmpersandError,fortComment,fortContComment,fortContError,fortContMark,fortDelimiter,fortExtended,fortFloat,fortHollerith,fortIntrinsic,fortNumber,fortOperator,fortParen,fortRightMargin,fortSpecific,fortString,fortTab,fortFmtString,fortIOControlSpec

" cpp is often used with Fortran
syn match cPreProc	"^\s*#\s*\(define\|elif\|else\|endif\|if\|ifdef\|\)\>.*"
syn match cPreProc	"^\s*#\s*\(ifndef\|include\|undef\)\>.*"

syn sync lines=20

if !exists("did_fortran_syntax_inits")
  let did_fortran_syntax_inits = 1
  " The default methods for highlighting.  Can be overridden later
  hi link fortTodo		Todo
  hi link fortTab		Todo
  hi link fortNumber		Number
  hi link fortHollerith		fortString
  hi link fortFloat		Float
  hi link fortOperator		Operator
  hi link fortLogical		Constant
  hi link fortStatement		Statement
  hi link fortRepeat		Repeat
  hi link fortGoto		fortRepeat
  hi link fortType		Type
  hi link fortDeclaration	fortType
  hi link fortUnitHeader	PreCondit
  hi link fortConditional	Conditional
  hi link fortIO		fortIntrinsic
  hi link fortIOControlSpec	fortIntrinsic
  hi link fortIntrinsic		Identifier
  hi link fortSpecific		fortIntrinsic
  hi link fortMemory		Statement
  hi link fortInclude		PreProc
  hi link fortExtended		Special
  hi link fortLabel		Special
  hi link fortLabelError	Error
  hi link fortComment		Comment
  hi link fortDocComment        DocComment
  hi link fortFirstLineMark	Conditional
  hi link fortContError		Error
  hi link fortContComment	fortComment
  hi link fortAmpersandError	Error
  hi link fortContMark		Todo
  hi link fortContMarkString	fortContMark
  hi link fortContMarkNext	fortContMark
  hi link fortStringInit1	fortString
  hi link fortStringInit2	fortString
  hi link fortString		String
  hi link fortFmtString		fortString
  hi link fortFormat		fortIO
  hi link fortParenError	Error
  hi link fortRightMargin	fortComment
  hi link cPreProc		PreProc

  " optional highlighting
  "hi link fortIdentifier	Identifier
  hi link fortDelimiter		PreProc
  "hi link fortSpecial		Special
endif

let b:current_syntax = "fortran90/95"

set nocindent
set smartindent

" ABB: print -> write (*,*)
iab pr\ write (*,*)

" ABB: module
iab module\ !! This module<CR>!! \author{Sixten Boeck}<CR>!! \date{today}<CR>module<CR><CR><TAB>use<CR><CR>implicit none<CR>! public types<CR>! public functions<CR>! public subroutines<CR>! interfaces<CR><CR>private<CR><CR><HOME>contains<CR><CR><HOME>end module<ESC>4GA

" ABB: function
iab function\ function() result ()<CR><TAB>use<CR>implicit none<CR><CR><BS>end function<CR><ESC>:-5<CR>wi

" ABB: subroutine
iab subroutine\ subroutine()<CR><TAB>use<CR>implicit none<CR><CR><BS>end subroutine<CR><ESC>:-5<CR>wi

" ABB: interface
iab interface\ interface<CR><TAB>module procedure<CR><BS>end interface<ESC>:-2<CR>A

" ABB: do
iab do\ do<CR><BS>end do<ESC>:-1<CR>A

" ABB: type
iab type\ type<CR><CR>end type<ESC>:-2<CR>A


" ABB: cvs
:iab cvs\ ! ------------------------------------------------------------------------<CR><HOME>!<CR><HOME>! $Id: $<CR><HOME>!<CR><HOME>! $Log: $<CR><HOME>!<CR><HOME>! ------------------------------------------------------------------------<CR><HOME><ESC>i


"EOF	vim: ts=8 noet tw=120 sw=8 sts=0
