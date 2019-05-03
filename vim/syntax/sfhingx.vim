" ---------------------------------------------------------------------------
" VIM syntax highlighting file for SFHIngX
" ---------------------------------------------------------------------------
" Installation:
"
"   (1) create file $HOME/.vim/syntax/syntax.vim containing
"       the following text:
"
"          augroup filetypedetect
"          au! BufRead,BufNewFile *.sx   setfiletype sfhingx
"          augroup END
"
"   (2) If you don't have a $HOME/.vim/syntax folder yet create it now
"
"         cd ~
"         mkdir -p .vim/syntax
"
"   (3) Copy or link this file to $HOME/.vim/syntax/sfhingx.vim
"
"         cd ~
"         cd .vim/syntax
"
"         cp     <YOUR_SFHINGX_PATH>/add-ons/sfhingx.vim .
"       OR
"         ln -sf <YOUR_SFHINGX_PATH>/add-ons/sfhingx.vim .
"
"   (4) Add to $HOME/.vimrc
"
"         so $HOME/.vim/syntax/syntax.vim
"         sy on
"
" ---------------------------------------------------------------------------
" Authors: Sixten Boeck         <boeck@sfhingx.de>
"          Christoph Freysoldt  <freyso@fhi-berlin.mpg.de>
" Date:    20/04/03
" ---------------------------------------------------------------------------

syn clear

" syn case ignore

syn keyword sxCommand     format
syn keyword sxCommand     sfhingx nextgroup=sxReqSemicolon
syn keyword sxCommand     include
syn keyword sxCommand     needs requires global

" --- structure
syn keyword sxGroup       structure species atom symmetry operator
syn keyword sxVariable    element cell
syn keyword sxVariable    operator 
syn keyword sxVariable    potential name element valenceCharge
syn keyword sxVariable    lMax lLoc lcaoOrbitals atomicRhoOcc
syn keyword sxVariable    rGauss reciprocalMass dampingMass ionicMass
syn keyword sxAttrib      movable relative nextgroup=sxReqSemicolon

" --- basis
syn keyword sxGroup       basis kPoint kPoints from to
syn keyword sxVariable    kUnits folding eCut mesh meshAccuracy chiEcut
syn keyword sxVariable    coords weight nPoints label
syn keyword sxAttrib      keepWavesOnDisc nextgroup=sxReqSemicolon

" --- Hamiltonian
syn keyword sxGroup       Hamiltonian 
syn keyword sxVariable    nEmptyStates nExcessElectrons ekt xc
syn keyword sxAttrib      LDA PBE PBE_LDA EXX READ_VXX nextgroup=sxReqSemicolon
syn keyword sxAttrib      spinPolarized nextgroup=sxReqSemicolon
" these are not followed necessarily by a semicolon, so no 
" nextgroup=sxReqSemicolon
syn keyword sxAttrib      CALC_ALL CALC_NONE CALC_KIN CALC_V_HARTREE
syn keyword sxAttrib      CALC_V_X CALC_V_C CALC_V_LOC CALC_UPDATE_RHO
syn keyword sxAttrib      CALC_V_NL CALC_V_XC CALC_V_EFF
syn keyword sxAttrib      EXX_WRITE_ALL EXX_WRITE_NONE EXX_WRITE_CHI_MATRIX
syn keyword sxAttrib      EXX_WRITE_CHI_ROWDIAG EXX_WRITE_E_G EXX_WRITE_VXR
syn keyword sxAttrib      EXX_WRITE_VXG

" --- initialGuess
syn keyword sxGroup       initialGuess waves rho lcao occ
syn keyword sxVariable    file maxSteps rhoMixing spinMoment
syn keyword sxAttrib      fromWaves random atomicOrbitals nextgroup=sxReqSemicolon

" --- main.elecMinim
syn keyword sxGroup       main
syn keyword sxGroup       SD WS DJ CCG DIIS_CCG DIAG
syn keyword sxVariable    maxSteps dEnergy dPsi printSteps
syn keyword sxVariable    deltaT gamma
syn keyword sxVariable    maxStepsCCG dEpsCCG
syn keyword sxVariable    nPulaySteps kerkerScaling kerkerDamping
syn keyword sxVariable    hContrib mixingMethod spinMixing
syn keyword sxAttrib      LINEAR PRECOND_LINEAR PULAY PRECOND_PULAY nextgroup=sxReqSemicolon
syn keyword sxAttrib      s p d f
syn keyword sxAttrib      keepRhoFixed keepOccFixed keepSpinFixed nextgroup=sxReqSemicolon
syn keyword sxAttrib      useFullBasis calcForces nextgroup=sxReqSemicolon
syn keyword sxVariable    eCutDiag nStatesDiag

syn keyword sxGroup       bandStructure subspaceDiag
syn keyword sxAttrib      autoSteps verbose printResidue
syn keyword sxVariable    maxSize overlap dEps sacrifyStates

" --- main.dampedNewton, quasiNewton, frozenPhonon, molDyn, synchronousTransit

syn keyword sxGroup      initStructure initHessian output 
syn keyword sxGroup      convergence dofRange performance
syn keyword sxGroup      initHistory randomVel devAtoms
syn keyword sxGroup      integrator thermostat
syn keyword sxGroup      initialStructure finalStructure  dofRange  


syn keyword sxVariable   initIdentity  file  samplePoints  saveStructure  
syn keyword sxVariable   saveWaves	saveHist    freezeRot   freezeIt 
syn keyword sxVariable   constraints     dAvgForceComponent  dEnergyStruct  
syn keyword sxVariable   dMaxForceComponent  maxStructSteps
syn keyword sxVariable   dAvgStructComponent   dMaxStructComponent 
syn keyword sxVariable   saveHessian  initEkin  gas deltaE  dof 
syn keyword sxVariable   temperature scheme  order  dt   timeSteps   
syn keyword sxVariable   deviation startDof  endDof  

syn keyword sxAttrib     pass optimizeRho expertOutput thirdOrderCor  nextgroup=sxReqSemicolon
 
syn keyword sxAttrib     writeHist   shiftToSticks	rescale   nextgroup=sxReqSemicolon
           
syn keyword sxAttrib     extrapolateWaves nextgroup=sxReqSemicolon
  
syn keyword sxGroup       dampedNewton quasiNewton molDyn frozenPhonon
syn keyword sxGroup       synchronousTransit 


" --- EXX specific keywords
syn keyword sxGroup     EXX_LOOP relaxRho
syn keyword sxVariable  writeControl
syn keyword sxAttrib    restart nextgroup=sxReqSemicolon


" --- isixServer and user
syn keyword sxGroup       isixServer user
syn keyword sxVariable    host port uid passwd


" --- numbers
syn match  sxNumber "[-+]\=\(\<\d[[:digit:]_]*L\=\>\|0[xX]\x[[:xdigit:]_]*\>\)"
syn match  sxNumber "[-+]\=\<\d[[:digit:]_]*[eE][\-+]\=\d\+"
syn match  sxNumber "[-+]\=\<\d[[:digit:]_]*\.[[:digit:]_]*\([eE][\-+]\=\d\+\)\="
syn match  sxNumber "[-+]\=\<\.[[:digit:]_]\+\([eE][\-+]\=\d\+\)\="

" --- some constants defined in parameters.sx
syn keyword sxNumber ANGSTROEM eV_by_Hartree EV_BY_HARTREE

" --- logicals
"syn match sxLogical "\(true\|TRUE\|false\|FALSE\|yes\|YES\|no\|NO\)"
syn keyword sxLogical true TRUE false FALSE yes YES no NO

" --- C++, Fortran, and shell like comments
syn region sxString start=/"/ end=/"/
syn region sxString start=/</ end=/>/ 
syn match  sxComment /#.*$/
syn match  sxComment /!.*$/
syn match  sxComment "//.*"


" --- C like comments
syn region sxComment    matchgroup=sxCommentStart start="/\*" matchgroup=NONE end="\*/" 
syntax match sxError    "\*/"

" --- check that each '=' is followed by "*;"
syn match sxError /=[^;[]*$/hs=e display
syn match sxError /=[^;[}]*}\{-1,}/hs=e display
syn match sxError "][^;,]*$" display

" --- check for missing semicolons ("=" is allowed, however)
syn match sxReqSemicolon "\_[ \n]*[^ ;=]" contained
hi link sxReqSemicolon sxError

" --- check brackets
"syn region sxBracket start="{" end="}" contains=ALLBUT,sxReqSemicolon,sxBrError
"syn match sxBrError "}" containedin=ALLBUT,sxBracket
"hi link sxBrError sxError


" ---- cross linking
hi link sxCommand       Include    
hi link sxGroup         Statement
hi link sxVariable      Type
hi link sxAttrib        SpecialChar
hi link sxChemElem      Identifier
hi link sxComment       Comment
hi link sxCommentStart  Comment
hi link sxNumber        Number
hi link sxLogical       Identifier
hi link sxString        String
hi link sxPath          String
hi link sxError         Error

let b:current_syntax = "sfhingx-input"
                                                 

" --- use always C indention for SFHIngX input files
set cindent															" CIndent...
set expandtab
set cinkeys=0{,0},:,0#,!,o,O,e,!<Tab>					" ...when...
set cinoptions="n2=10+17(0)100*100							" ...and how
