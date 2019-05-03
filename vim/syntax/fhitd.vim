" Syntax highliting file for FHI/ng - the next generation
"
" Author: Sixten Boeck <boeck@fhi-berlin.mpg.de>
" Date:   01/02/99

syn clear

" syn case ignore

syn keyword fhiReqGrp      slab cell
syn keyword fhiOptGrp      species atoms
syn keyword fhiReqGrp      coords tau
syn keyword fhiOptGrp      relax
syn keyword fhiAttribute   type nItems dim dims min max optional
syn keyword fhiAttribute   minItems maxItems nItems val

" syn match   fhiOptGrp     "perturb_.*$"

syn match fhiComment "#.*$"
"syn match fhiNumber  "\<\d\+\(u\=l\=\|lu\|f\)\>"
syn match fhiLogical	"\.\(true\|t\|false\|f\)\."


syn region fhiString start=/"/ end =/"/


" ---- cross linking

hi link fhiReqGrp        Statement
hi link fhiOptGrp        Identifier
hi link fhiAttribute     Identifier
hi link fhiChemElem      Type
hi link fhiComment       Comment
hi link fhiNumber        Number
hi link fhiLogical       Conditional
hi link fhiString        String

let b:current_syntax = "fhi-input"
						 

