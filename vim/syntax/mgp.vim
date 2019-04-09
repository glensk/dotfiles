syn clear

syn case ignore

syn match   mgpPercent   "^%"
syn match   mgpComment   "^%%.*$"
syn match   mgpNumber    "\<\d\+\(u\=l\=\|lu\|f\)\>"

syn keyword mgpStatement tab default page bgrad center left right 
syn keyword mgpStatement fore back size rcutin lcutin shrink pause
syn keyword mgpStatement hgap vgap fontleftfill bar vfont prefix
syn keyword mgpStatement image system
syn keyword mgpColor     white yellow green blue red cyan brown black


hi link     mgpComment   Comment
hi link     mgpStatement Statement
hi link     mgpPercent   Statement
hi link     mgpNumber    Number


