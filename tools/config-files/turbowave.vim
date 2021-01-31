" Vim syntax file
" Language:	turboWAVE
" Maintainer: Daniel Gordon

" Syntax highlighting for turboWAVE input files:
" Put turbowave.vim into ~/.vim/syntax/ and put filetype.vim into ~/.vim
" This causes stdin, stdin.txt, and *.tw to be highlighted
" To highlight any other file from within vim, use ``:set syntax=turbowave``

" Quit when a (custom) syntax file was already loaded
if exists("b:current_syntax")
  finish
endif

syn match twPreproc "#define"
syn match twPreproc "#include"
syn match twPreproc "#ifdef"
syn match twPreproc "#ifndef"
syn match twPreproc "#else"
syn match twPreproc "#endif"

syn keyword	twStatement	new generate get for

syn match twUnit "\v(\s|\(|\{|\,|\=|$)[\+\-]=([0-9]+\.=[0-9]*|\.[0-9]+)([eE][\+\-]=[0-9]+)=\s*\[(deg|rad|mrad|urad|cm2|m2|cm2/s|m2/s|um|mm|cm|m|fs|ps|ns|us|s|/m3|/cm3|J/m3|J/cm3|eV|K|V|webers/m|G\*cm|V/m|V/cm|T|G)\](\s|\)|\}|\,|\n)"hs=s+1,he=e-1,me=e-1
syn region twUserMacro display start="\$" end="\v(\s|$)"he=s-1,re=s-1,me=s-1
syn region twComment start="/\*" end="\*/"
syn region twComment start="//" end="$" keepend

syn keyword twconst true false on off yes no
syn keyword twconst none full
syn keyword twconst deterministic statistical variable fixed triggered maintained
syn keyword twconst absorbing emitting reflecting axisymmetric periodic
syn keyword twconst adk ppt mpi kyh pmpb
syn keyword twconst quintic sin2 sech
syn match twmode "hermite gauss"
syn match twmode "laguerre gauss"
syn match twmode "plane wave"
syn match twmode "airy disc"
syn match twmode "bessel beam"
syn match twmode "multipole"

syn region twString start="'" end="'"
syn region twString start='"' end='"'

syn region twBlock start="{" end="}" fold transparent contains=ALL

hi def link twPreproc Define
hi def link twStatement Statement
hi def link twUnit Macro
hi def link twUserMacro Macro
hi def link twComment Comment
hi def link twconst Constant
hi def link twmode Constant
hi def link twString String

let b:current_syntax = "turbowave"
